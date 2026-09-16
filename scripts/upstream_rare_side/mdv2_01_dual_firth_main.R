#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(args[[i + 1]])
  default
}
config_path <- get_arg("--config", "config/mdv2_dual_firth_gate_v1.yaml")
threads_arg <- suppressWarnings(as.integer(get_arg("--threads", NA_character_)))

required_pkgs <- c("yaml", "data.table", "brglm2", "parallel")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Missing R package(s): ", paste(missing_pkgs, collapse = ", "))

cfg <- yaml::read_yaml(config_path)
root <- normalizePath(cfg$project_root, mustWork = TRUE)
abs_path <- function(x) normalizePath(file.path(root, x), mustWork = FALSE)
threads <- if (is.finite(threads_arg) && threads_arg >= 1) threads_arg else as.integer(cfg$runtime$threads)
threads <- max(1L, threads)

final_pass <- abs_path(cfg$outputs$final_pass)
if (file.exists(final_pass)) {
  fp <- readLines(final_pass, warn = FALSE)
  if (!any(grepl("^status=PASS$", fp))) stop("Existing MDV2_PASS is malformed")
  message("[PASS] MDV2_PASS already exists; immutable freeze prevents re-analysis.")
  quit(save = "no", status = 0)
}
preflight_pass <- abs_path(cfg$outputs$preflight_pass)
if (!file.exists(preflight_pass)) stop("MDV2 preflight PASS is missing")
pf_lines <- readLines(preflight_pass, warn = FALSE)
if (!any(grepl("^status=PASS$", pf_lines))) stop("MDV2 preflight PASS is malformed")
pf_ver <- sub("^version=", "", pf_lines[grepl("^version=", pf_lines)])
if (length(pf_ver) != 1L || !identical(pf_ver, as.character(cfg$version))) {
  stop("MDV2 preflight version mismatch: observed=", paste(pf_ver, collapse = ","),
       "; expected=", cfg$version, ". Re-run mdv2_preflight with the current v1.1 contract before analysis.")
}

u <- data.table::fread(abs_path(cfg$inputs$gene_universe), sep = "\t", showProgress = FALSE)
pm <- data.table::fread(abs_path(cfg$inputs$pathway_membership), sep = "\t", showProgress = FALSE)
cc <- cfg$columns
alpha <- as.numeric(cfg$primary_model$alpha)
ci_level <- as.numeric(cfg$primary_model$ci_level)
zcrit <- stats::qnorm(1 - (1 - ci_level) / 2)

# Minimal internal assertions: detailed identity/integrity checks were already done in preflight.
stopifnot(nrow(u) == as.integer(cfg$expected_frozen_state$protein_coding_universe_n))
stopifnot(sum(u[[cc$group]] == "Overlap_raw") == as.integer(cfg$expected_frozen_state$R0_Overlap_raw_n))
stopifnot(sum(as.integer(u[[cc$coreseed]]) == 1L) == as.integer(cfg$expected_frozen_state$CoreSeed_n))

# Canonical internal pathway columns, then freeze one row per pathway. Keep the full 10-500 family,
# including pathways with no contrast-level predictor variation.
pmv <- pm[, .(
  source = get(cc$pathway_source), pathway_id = get(cc$pathway_id),
  pathway_name = get(cc$pathway_name), approved_symbol = get(cc$pathway_symbol)
)]
paths <- pmv[, .(
  pathway_name = unique(pathway_name)[[1]],
  global_pathway_size = data.table::uniqueN(approved_symbol)
), by = .(source, pathway_id)]
data.table::setorder(paths, source, pathway_id)

pm_key <- paste(pmv$source, pmv$pathway_id, sep = "\r")
member_sets <- split(pmv$approved_symbol, pm_key)
path_keys <- paste(paths$source, paths$pathway_id, sep = "\r")
stopifnot(all(path_keys %in% names(member_sets)))

normalize_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)

build_contrast <- function(contrast, mode = "primary") {
  group <- u[[cc$group]]
  core <- as.integer(u[[cc$coreseed]]) == 1L
  keep <- rep(FALSE, nrow(u))
  y <- rep(NA_integer_, nrow(u))

  if (mode == "coreseed") {
    if (contrast == "C1") {
      keep <- core | group == "IEI_only"
      y[keep] <- as.integer(core[keep])
    } else {
      keep <- core | group == "ASD_only"
      y[keep] <- as.integer(core[keep])
    }
  } else {
    if (contrast == "C1") {
      keep <- group %in% c("Overlap_raw", "IEI_only")
      y[keep] <- as.integer(group[keep] == "Overlap_raw")
    } else {
      keep <- group %in% c("Overlap_raw", "ASD_only")
      y[keep] <- as.integer(group[keep] == "Overlap_raw")
    }
  }

  if (mode == "mhc_exclusion") {
    chr <- normalize_chr(u[[cc$chr]])
    st <- as.numeric(u[[cc$start]]); en <- as.numeric(u[[cc$end]])
    in_mhc <- chr == as.character(cfg$mhc_grch37$chr) & en >= as.numeric(cfg$mhc_grch37$start) & st <= as.numeric(cfg$mhc_grch37$end)
    keep <- keep & !in_mhc
  }

  d <- u[keep]
  d[, outcome := y[keep]]
  d
}

formula_for <- function(length_mode = c("spline", "linear"), extra = NULL) {
  length_mode <- match.arg(length_mode)
  len_col <- cc$log_gene_length; gc_col <- cc$gc_z; ann_col <- cc$annotation_degree
  len_term <- if (length_mode == "spline") sprintf("splines::ns(%s, df=%d)", len_col, as.integer(cfg$primary_model$spline_df)) else len_col
  rhs <- c("pathway_member", len_term, gc_col, sprintf("scale(log1p(%s))", ann_col))
  if (!is.null(extra)) rhs <- c(rhs, extra)
  stats::as.formula(paste("outcome ~", paste(rhs, collapse = " + ")))
}

fit_one <- function(d, member_symbols, formula, contrast, source, pathway_id, pathway_name, global_size) {
  d <- data.table::copy(d)
  d[, pathway_member := as.integer(get(cc$universe_symbol) %in% member_symbols)]
  n_case <- sum(d$outcome == 1L)
  n_control <- sum(d$outcome == 0L)
  pm_case <- sum(d$pathway_member[d$outcome == 1L])
  pm_control <- sum(d$pathway_member[d$outcome == 0L])

  base <- list(
    source = source, pathway_id = pathway_id, pathway_name = pathway_name,
    global_pathway_size = as.integer(global_size), contrast = contrast,
    n_case = as.integer(n_case), n_control = as.integer(n_control),
    pathway_case = as.integer(pm_case), pathway_control = as.integer(pm_control)
  )

  if (length(unique(d$pathway_member)) < 2L) {
    return(c(base, list(beta = NA_real_, SE = NA_real_, OR = NA_real_, CI95L = NA_real_, CI95U = NA_real_, P = 1.0,
                        fit_status = "NO_VARIATION", warning = "pathway_member has no variation in contrast subset")))
  }

  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm(formula, data = d, family = stats::binomial("logit"), method = brglm2::brglmFit,
                 type = cfg$primary_model$brglm2_type),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(c(base, list(beta = NA_real_, SE = NA_real_, OR = NA_real_, CI95L = NA_real_, CI95U = NA_real_, P = 1.0,
                        fit_status = "FIT_ERROR", warning = conditionMessage(fit))))
  }

  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
  vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  if (is.null(cf) || is.null(vc) || !("pathway_member" %in% names(cf)) || !("pathway_member" %in% rownames(vc))) {
    return(c(base, list(beta = NA_real_, SE = NA_real_, OR = NA_real_, CI95L = NA_real_, CI95U = NA_real_, P = 1.0,
                        fit_status = "COEF_MISSING", warning = paste(c(warnings, "pathway coefficient/vcov missing"), collapse = " | "))))
  }
  beta <- as.numeric(cf[["pathway_member"]])
  se <- sqrt(as.numeric(vc["pathway_member", "pathway_member"]))
  p <- if (is.finite(beta) && is.finite(se) && se > 0) 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE) else NA_real_
  if (!is.finite(beta) || !is.finite(se) || se <= 0 || !is.finite(p) || p < 0 || p > 1) {
    return(c(base, list(beta = beta, SE = se, OR = NA_real_, CI95L = NA_real_, CI95U = NA_real_, P = 1.0,
                        fit_status = "COEF_INVALID", warning = paste(c(warnings, "non-finite/invalid coefficient statistics"), collapse = " | "))))
  }
  c(base, list(beta = beta, SE = se, OR = exp(beta), CI95L = exp(beta - zcrit * se), CI95U = exp(beta + zcrit * se), P = p,
               fit_status = "PASS", warning = paste(unique(warnings), collapse = " | ")))
}

run_dual <- function(analysis_id, mode = "primary", length_mode = "spline", extra = NULL) {
  d1 <- build_contrast("C1", mode = mode)
  d2 <- build_contrast("C2", mode = mode)
  f <- formula_for(length_mode = length_mode, extra = extra)

  worker <- function(i) {
    key <- path_keys[[i]]
    ms <- member_sets[[key]]
    p <- paths[i]
    r1 <- fit_one(d1, ms, f, "C1", p$source, p$pathway_id, p$pathway_name, p$global_pathway_size)
    r2 <- fit_one(d2, ms, f, "C2", p$source, p$pathway_id, p$pathway_name, p$global_pathway_size)
    list(r1 = r1, r2 = r2)
  }
  idx <- seq_len(nrow(paths))
  rr <- if (.Platform$OS.type == "unix" && threads > 1L) parallel::mclapply(idx, worker, mc.cores = threads, mc.preschedule = TRUE) else lapply(idx, worker)
  c1 <- data.table::rbindlist(lapply(rr, `[[`, "r1"), fill = TRUE)
  c2 <- data.table::rbindlist(lapply(rr, `[[`, "r2"), fill = TRUE)

  stopifnot(nrow(c1) == nrow(paths), nrow(c2) == nrow(paths))
  out <- data.table::data.table(
    analysis_id = analysis_id,
    source = paths$source, pathway_id = paths$pathway_id, pathway_name = paths$pathway_name,
    global_pathway_size = paths$global_pathway_size,
    C1_n_case = c1$n_case, C1_n_control = c1$n_control, C1_pathway_case = c1$pathway_case, C1_pathway_control = c1$pathway_control,
    beta_C1 = c1$beta, SE_C1 = c1$SE, OR_C1 = c1$OR, CI95L_C1 = c1$CI95L, CI95U_C1 = c1$CI95U, P_C1 = c1$P,
    fit_status_C1 = c1$fit_status, warning_C1 = c1$warning,
    C2_n_case = c2$n_case, C2_n_control = c2$n_control, C2_pathway_case = c2$pathway_case, C2_pathway_control = c2$pathway_control,
    beta_C2 = c2$beta, SE_C2 = c2$SE, OR_C2 = c2$OR, CI95L_C2 = c2$CI95L, CI95U_C2 = c2$CI95U, P_C2 = c2$P,
    fit_status_C2 = c2$fit_status, warning_C2 = c2$warning
  )
  out[, P_conj := pmax(P_C1, P_C2)]
  out[, q_conj := stats::p.adjust(P_conj, method = cfg$primary_model$fdr_method), by = source]
  out[, direction_pass := is.finite(OR_C1) & is.finite(OR_C2) & OR_C1 > 1 & OR_C2 > 1]
  out[, conjunction_testable := fit_status_C1 != "NO_VARIATION" & fit_status_C2 != "NO_VARIATION"]
  out[, primary_pass := fit_status_C1 == "PASS" & fit_status_C2 == "PASS" & direction_pass & q_conj < alpha]
  out[]
}

message("[INFO] Running primary dual bias-reduced logistic conjunction across full frozen pathway families...")
primary <- run_dual("PRIMARY_RAW_OVERLAP", mode = "primary", length_mode = "spline")

# Required sensitivities.
sens_list <- list()
registry <- data.table::data.table(analysis_id = character(), required = logical(), status = character(), detail = character())
for (s in cfg$sensitivities$required) {
  sid <- s$id; mode <- s$mode
  message("[INFO] Running required sensitivity: ", sid)
  if (mode == "coreseed") {
    z <- run_dual(sid, mode = "coreseed", length_mode = "spline")
  } else if (mode == "linear_length") {
    z <- run_dual(sid, mode = "primary", length_mode = "linear")
  } else if (mode == "mhc_exclusion") {
    z <- run_dual(sid, mode = "mhc_exclusion", length_mode = "spline")
  } else {
    stop("Unknown required sensitivity mode: ", mode)
  }
  sens_list[[sid]] <- z
  registry <- data.table::rbindlist(list(registry, data.table::data.table(analysis_id = sid, required = TRUE, status = "COMPLETED", detail = s$description)), use.names = TRUE)
}

# Optional structural covariate sensitivities are run only if the frozen universe actually contains the configured covariate.
if (length(cfg$sensitivities$optional_covariates)) {
  for (s in cfg$sensitivities$optional_covariates) {
    sid <- s$id; col <- s$column; tr <- s$transform
    if (!(col %in% names(u))) {
      registry <- data.table::rbindlist(list(registry, data.table::data.table(analysis_id = sid, required = isTRUE(s$required), status = "SKIPPED_NOT_IN_FROZEN_MDV1", detail = col)), use.names = TRUE)
      next
    }
    xv <- suppressWarnings(as.numeric(u[[col]]))
    if (length(xv) != nrow(u) || any(!is.finite(xv)) || any(xv < 0 & identical(tr, "scale(log1p)"))) {
      registry <- data.table::rbindlist(list(registry, data.table::data.table(analysis_id = sid, required = isTRUE(s$required), status = "SKIPPED_INVALID_FROZEN_COLUMN", detail = col)), use.names = TRUE)
      next
    }
    extra <- if (identical(tr, "scale(log1p)")) sprintf("scale(log1p(%s))", col) else if (identical(tr, "scale")) sprintf("scale(%s)", col) else stop("Unsupported optional transform: ", tr)
    message("[INFO] Running optional covariate sensitivity: ", sid)
    sens_list[[sid]] <- run_dual(sid, mode = "primary", length_mode = "spline", extra = extra)
    registry <- data.table::rbindlist(list(registry, data.table::data.table(analysis_id = sid, required = isTRUE(s$required), status = "COMPLETED", detail = paste(col, tr))), use.names = TRUE)
  }
}
sens <- data.table::rbindlist(sens_list, use.names = TRUE, fill = TRUE)

# Model QC registry. NO_VARIATION is structural and not a fit failure; any other non-PASS status is a technical error.
# v1.1 HOTFIX: never place a length-1 constant inside data.table `by=`. Recent data.table
# versions require every explicit `by` item to have nrow(dt) length. Aggregate first, then
# assign the scalar contrast label. This is semantically identical and version-robust.
qc_one_contrast <- function(dt, status_col, contrast_label) {
  required <- c("analysis_id", "source", status_col)
  miss <- setdiff(required, names(dt))
  if (length(miss)) stop("QC aggregation missing column(s): ", paste(miss, collapse = ", "))
  z <- dt[, .N, by = c("analysis_id", "source", status_col)]
  data.table::setnames(z, status_col, "fit_status")
  z[, contrast := contrast_label]
  data.table::setcolorder(z, c("analysis_id", "source", "contrast", "fit_status", "N"))
  z[]
}
qc_one <- function(dt) {
  data.table::rbindlist(list(
    qc_one_contrast(dt, "fit_status_C1", "C1"),
    qc_one_contrast(dt, "fit_status_C2", "C2")
  ), use.names = TRUE, fill = TRUE)
}
qc <- data.table::rbindlist(list(qc_one(primary), qc_one(sens)), use.names = TRUE, fill = TRUE)
if (!nrow(qc)) stop("QC aggregation unexpectedly returned zero rows")
if (anyNA(qc$analysis_id) || anyNA(qc$source) || anyNA(qc$contrast) || anyNA(qc$fit_status) || anyNA(qc$N)) {
  stop("QC aggregation produced missing key/count values")
}
data.table::setorder(qc, analysis_id, source, contrast, fit_status)

write_tsv_atomic <- function(dt, path) {
  tmp <- paste0(path, ".tmp")
  data.table::fwrite(dt, tmp, sep = "\t", quote = FALSE, na = "NA")
  if (!file.rename(tmp, path)) stop("Failed to atomically promote: ", path)
}
write_tsvgz_atomic <- function(dt, path) {
  tmp <- paste0(path, ".tmp")
  con <- gzfile(tmp, open = "wt")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  utils::write.table(dt, con, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = "NA")
  close(con)
  if (!file.rename(tmp, path)) stop("Failed to atomically promote: ", path)
}

primary_path <- abs_path(cfg$outputs$primary_results)
hits_path <- abs_path(cfg$outputs$primary_hits)
sens_path <- abs_path(cfg$outputs$sensitivity_results)
registry_path <- abs_path(cfg$outputs$sensitivity_registry)
qc_path <- abs_path(cfg$outputs$model_qc)
complete_path <- abs_path(cfg$outputs$analysis_complete)

dir.create(dirname(primary_path), recursive = TRUE, showWarnings = FALSE)
write_tsvgz_atomic(primary, primary_path)
write_tsv_atomic(primary[primary_pass == TRUE], hits_path)
write_tsvgz_atomic(sens, sens_path)
write_tsv_atomic(registry, registry_path)
write_tsv_atomic(qc, qc_path)

n_hits <- sum(primary$primary_pass)
complete <- c(
  "MDV2_ANALYSIS_COMPLETE",
  "status=COMPLETE",
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("version=", cfg$version),
  paste0("threads=", threads),
  paste0("primary_pathways_n=", nrow(primary)),
  paste0("primary_hits_n=", n_hits),
  paste0("required_sensitivities_n=", sum(registry$required & registry$status == "COMPLETED")),
  paste0("optional_sensitivities_completed_n=", sum(!registry$required & registry$status == "COMPLETED")),
  paste0("optional_sensitivities_skipped_n=", sum(!registry$required & grepl("^SKIPPED", registry$status))),
  paste0("R_version=", R.version.string),
  paste0("brglm2_version=", as.character(utils::packageVersion("brglm2"))),
  paste0("data.table_version=", as.character(utils::packageVersion("data.table"))),
  "GWAS_inputs_used=NO",
  "next=MDV2_POSTFLIGHT_VALIDATION"
)
writeLines(complete, paste0(complete_path, ".tmp"))
if (!file.rename(paste0(complete_path, ".tmp"), complete_path)) stop("Failed to promote analysis complete sentinel")
message("[PASS] MDV-2 primary + sensitivity analysis completed. Postflight gate must still pass before MDV-2 closes.")
