suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(yaml)
})

mdv4_stop <- function(...) stop(..., call. = FALSE)
mdv4_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
mdv4_dir <- function(path) { dir.create(path, recursive = TRUE, showWarnings = FALSE); invisible(path) }

mdv4_read_cfg <- function(path) {
  if (!file.exists(path)) mdv4_stop("Missing config: ", path)
  yaml::read_yaml(path)
}

mdv4_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
        out[[key]] <- TRUE
        i <- i + 1L
      } else {
        out[[key]] <- args[[i + 1L]]
        i <- i + 2L
      }
    } else i <- i + 1L
  }
  out
}

mdv4_atomic_fwrite <- function(x, path, ...) {
  mdv4_dir(dirname(path))
  tmp <- paste0(path, ".tmp")
  data.table::fwrite(x, tmp, sep = "\t", quote = FALSE, ...)
  if (!file.rename(tmp, path)) mdv4_stop("Atomic rename failed for ", path)
  invisible(path)
}

mdv4_atomic_lines <- function(lines, path) {
  mdv4_dir(dirname(path))
  tmp <- paste0(path, ".tmp")
  writeLines(as.character(lines), tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) mdv4_stop("Atomic rename failed for ", path)
  invisible(path)
}

mdv4_kv <- function(path) {
  if (!file.exists(path)) mdv4_stop("Missing sentinel: ", path)
  z <- readLines(path, warn = FALSE)
  z <- z[grepl("=", z, fixed = TRUE)]
  k <- sub("=.*$", "", z)
  v <- sub("^[^=]*=", "", z)
  stats::setNames(v, k)
}

mdv4_sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  z <- suppressWarnings(system2("sha256sum", path, stdout = TRUE, stderr = TRUE))
  if (length(z) < 1L) return(NA_character_)
  strsplit(z[[1]], "[[:space:]]+")[[1]][1]
}

mdv4_verify_checksum_manifest <- function(manifest, root = ".") {
  if (!file.exists(manifest)) return(list(ok = FALSE, output = paste("missing", manifest)))
  old <- getwd(); on.exit(setwd(old), add = TRUE); setwd(root)
  z <- suppressWarnings(system2("sha256sum", c("-c", manifest), stdout = TRUE, stderr = TRUE))
  status <- attr(z, "status"); if (is.null(status)) status <- 0L
  list(ok = identical(as.integer(status), 0L), output = paste(z, collapse = "\n"))
}

mdv4_length_quintile <- function(log_length) {
  r <- data.table::frank(as.numeric(log_length), ties.method = "average")
  pmin(5L, pmax(1L, floor((r - 1) / length(log_length) * 5) + 1L))
}

mdv4_prepare_universe <- function(u, cfg) {
  cc <- cfg$columns
  req <- unlist(cc[c("hgnc_id", "hgnc_symbol", "coreseed", "log_gene_length", "gc_z", "annotation_z")])
  miss <- setdiff(req, names(u))
  if (length(miss)) mdv4_stop("Universe missing columns: ", paste(miss, collapse = ", "))
  for (nm in c(cc$log_gene_length, cc$gc_z, cc$annotation_z)) {
    if (any(!is.finite(as.numeric(u[[nm]])))) mdv4_stop("Non-finite values in ", nm)
  }
  x <- as.numeric(u[[cc$log_gene_length]])
  u[, match_loglen_z := as.numeric(scale(x))]
  u[, match_gc_z := as.numeric(get(cc$gc_z))]
  u[, match_annotation_z := as.numeric(get(cc$annotation_z))]
  u[, length_quintile := mdv4_length_quintile(get(cc$log_gene_length))]
  u
}

mdv4_balance_smd <- function(core_dt, pseudo_dt) {
  vars <- c("match_loglen_z", "match_gc_z", "match_annotation_z")
  out <- setNames(rep(NA_real_, length(vars)), vars)
  for (v in vars) {
    a <- as.numeric(core_dt[[v]]); b <- as.numeric(pseudo_dt[[v]])
    sp <- sqrt(((length(a)-1)*stats::var(a) + (length(b)-1)*stats::var(b)) / (length(a)+length(b)-2))
    out[[v]] <- if (is.finite(sp) && sp > 0) (mean(a) - mean(b)) / sp else NA_real_
  }
  out
}

mdv4_build_candidate_pools <- function(u, cfg, pool_size = NULL) {
  cc <- cfg$columns
  if (is.null(pool_size)) pool_size <- as.integer(cfg$matching$pool_size_primary)
  core <- u[get(cc$coreseed) == 1L]
  noncore <- u[get(cc$coreseed) == 0L]
  if (nrow(core) != as.integer(cfg$expected$coreseed_n)) mdv4_stop("CoreSeed count mismatch")
  res <- vector("list", nrow(core))
  for (ii in seq_len(nrow(core))) {
    s <- core[ii]
    q0 <- as.integer(s$length_quintile)
    cand <- copy(noncore)
    cand[, qdist := abs(length_quintile - q0)]
    cand[, distance := sqrt(
      (match_loglen_z - s$match_loglen_z)^2 +
      (match_gc_z - s$match_gc_z)^2 +
      (match_annotation_z - s$match_annotation_z)^2
    )]
    # Primary priority: same quintile; fallback: adjacent quintiles by quintile distance, then multivariate distance.
    setorderv(cand, c("qdist", "distance", cc$hgnc_symbol))
    cand <- cand[seq_len(min(pool_size, nrow(cand)))]
    if (nrow(cand) < pool_size) mdv4_stop("Insufficient candidates for seed ", s[[cc$hgnc_symbol]])
    res[[ii]] <- data.table(
      seed_HGNC_id = s[[cc$hgnc_id]],
      seed_HGNC_symbol = s[[cc$hgnc_symbol]],
      seed_length_quintile = q0,
      seed_match_loglen_z = s$match_loglen_z,
      seed_match_gc_z = s$match_gc_z,
      seed_match_annotation_z = s$match_annotation_z,
      candidate_HGNC_id = cand[[cc$hgnc_id]],
      candidate_HGNC_symbol = cand[[cc$hgnc_symbol]],
      candidate_length_quintile = cand$length_quintile,
      candidate_match_loglen_z = cand$match_loglen_z,
      candidate_match_gc_z = cand$match_gc_z,
      candidate_match_annotation_z = cand$match_annotation_z,
      candidate_rank = seq_len(nrow(cand)),
      distance = cand$distance,
      expansion_level = cand$qdist
    )
  }
  rbindlist(res)
}

mdv4_generate_assignments <- function(pools, B, seed, max_attempts = 1000L) {
  B <- as.integer(B); seed <- as.integer(seed); max_attempts <- as.integer(max_attempts)
  seeds <- unique(pools$seed_HGNC_symbol)
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  pool_list <- split(pools$candidate_HGNC_symbol, pools$seed_HGNC_symbol)
  accepted_keys <- new.env(parent = emptyenv(), hash = TRUE)
  out <- vector("list", B)
  redraw_total <- 0L
  for (b in seq_len(B)) {
    success <- FALSE
    for (attempt in seq_len(max_attempts)) {
      order_seed <- sample(seeds, length(seeds), replace = FALSE)
      chosen <- setNames(rep(NA_character_, length(seeds)), seeds)
      used <- character()
      dead <- FALSE
      for (s in order_seed) {
        avail <- setdiff(pool_list[[s]], used)
        if (!length(avail)) { dead <- TRUE; break }
        g <- sample(avail, 1L)
        chosen[[s]] <- g
        used <- c(used, g)
      }
      if (dead) { redraw_total <- redraw_total + 1L; next }
      key <- paste(sort(unname(chosen)), collapse = "|")
      if (exists(key, envir = accepted_keys, inherits = FALSE)) {
        redraw_total <- redraw_total + 1L; next
      }
      assign(key, TRUE, envir = accepted_keys)
      out[[b]] <- data.table(set_id = b, seed_HGNC_symbol = seeds, matched_HGNC_symbol = unname(chosen[seeds]))
      success <- TRUE
      break
    }
    if (!success) mdv4_stop("Failed to generate unique pseudo-set ", b, " after ", max_attempts, " attempts")
  }
  ans <- rbindlist(out)
  attr(ans, "redraw_total") <- redraw_total
  ans
}

mdv4_assignments_wide <- function(assignments) {
  dcast(assignments, set_id ~ seed_HGNC_symbol, value.var = "matched_HGNC_symbol")
}

mdv4_assignment_balance <- function(assignments, u, cfg) {
  cc <- cfg$columns
  core <- u[get(cc$coreseed) == 1L]
  covmap <- u[, .(HGNC_symbol = get(cc$hgnc_symbol), match_loglen_z, match_gc_z, match_annotation_z, length_quintile)]
  setkey(covmap, HGNC_symbol)
  ids <- sort(unique(assignments$set_id))
  out <- vector("list", length(ids))
  core_q <- table(factor(core$length_quintile, levels = 1:5))
  for (ii in seq_along(ids)) {
    sid <- ids[[ii]]
    genes <- assignments[set_id == sid, matched_HGNC_symbol]
    p <- covmap[J(genes)]
    if (nrow(p) != length(genes) || anyNA(p$HGNC_symbol)) mdv4_stop("Balance join failed for set ", sid)
    smd <- mdv4_balance_smd(core, p)
    pq <- table(factor(p$length_quintile, levels = 1:5))
    out[[ii]] <- data.table(
      set_id = sid,
      n = length(genes),
      unique_n = uniqueN(genes),
      smd_log_gene_length = smd[["match_loglen_z"]],
      smd_GC = smd[["match_gc_z"]],
      smd_annotation = smd[["match_annotation_z"]],
      max_abs_smd = max(abs(smd), na.rm = TRUE),
      exact_length_quintile_profile = identical(as.integer(pq), as.integer(core_q))
    )
  }
  rbindlist(out)
}

mdv4_logor <- function(a, k, m, N, correction = 0.5) {
  a <- as.numeric(a); k <- as.numeric(k); m <- as.numeric(m); N <- as.numeric(N)
  b <- k - a
  c <- m - a
  d <- N - k - c
  z <- (a == 0 | b == 0 | c == 0 | d == 0)
  aa <- a + ifelse(z, correction, 0)
  bb <- b + ifelse(z, correction, 0)
  cc <- c + ifelse(z, correction, 0)
  dd <- d + ifelse(z, correction, 0)
  log((aa * dd) / (bb * cc))
}

mdv4_make_membership_matrix <- function(u, mem, paths, cfg) {
  cc <- cfg$columns
  pkey <- paste(paths$source, paths$pathway_id, sep = "||")
  mem[, path_key := paste(source, pathway_id, sep = "||")]
  ii <- match(mem[[cfg$membership$hgnc_id]], u[[cc$hgnc_id]])
  jj <- match(mem$path_key, pkey)
  keep <- !is.na(ii) & !is.na(jj)
  M <- Matrix::sparseMatrix(i = ii[keep], j = jj[keep], x = 1,
                            dims = c(nrow(u), nrow(paths)),
                            dimnames = list(u[[cc$hgnc_symbol]], pkey))
  M
}

mdv4_empirical_from_assignments <- function(assignments, u, M, paths, cfg, B_expected = NULL, chunk_size = NULL) {
  cc <- cfg$columns
  ids <- sort(unique(assignments$set_id))
  if (!is.null(B_expected) && length(ids) != as.integer(B_expected)) mdv4_stop("Pseudo-set B mismatch")
  if (is.null(chunk_size)) chunk_size <- as.integer(cfg$empirical$chunk_size)
  N <- nrow(u); k <- as.integer(cfg$expected$coreseed_n)
  core_idx <- which(u[[cc$coreseed]] == 1L)
  m <- as.numeric(Matrix::colSums(M))
  aobs <- as.numeric(Matrix::colSums(M[core_idx, , drop = FALSE]))
  Tobs <- mdv4_logor(aobs, k = k, m = m, N = N, correction = cfg$empirical$haldane_correction)
  exceed <- numeric(ncol(M)); sumT <- numeric(ncol(M)); sumT2 <- numeric(ncol(M))
  gene_index <- setNames(seq_len(nrow(u)), u[[cc$hgnc_symbol]])
  id_chunks <- split(ids, ceiling(seq_along(ids) / chunk_size))
  for (chunk_ids in id_chunks) {
    aa <- assignments[set_id %in% chunk_ids]
    ri <- match(aa$set_id, chunk_ids)
    cj <- unname(gene_index[aa$matched_HGNC_symbol])
    if (anyNA(cj)) mdv4_stop("Unknown matched gene in assignments")
    S <- Matrix::sparseMatrix(i = ri, j = cj, x = 1,
                              dims = c(length(chunk_ids), nrow(u)))
    A <- as.matrix(S %*% M)
    MM <- matrix(m, nrow = nrow(A), ncol = length(m), byrow = TRUE)
    Bc <- k - A; Cc <- MM - A; Dc <- N - k - Cc
    z <- (A == 0 | Bc == 0 | Cc == 0 | Dc == 0)
    corr <- as.numeric(cfg$empirical$haldane_correction)
    Tm <- log(((A + corr*z) * (Dc + corr*z)) / ((Bc + corr*z) * (Cc + corr*z)))
    TO <- matrix(Tobs, nrow = nrow(Tm), ncol = length(Tobs), byrow = TRUE)
    exceed <- exceed + colSums(Tm >= (TO - as.numeric(cfg$empirical$comparison_tolerance)))
    sumT <- sumT + colSums(Tm)
    sumT2 <- sumT2 + colSums(Tm^2)
  }
  Bn <- length(ids)
  pemp <- (1 + exceed) / (Bn + 1)
  null_mean <- sumT / Bn
  null_var <- pmax(0, (sumT2 - Bn * null_mean^2) / pmax(1, Bn - 1))
  data.table(
    source = paths$source,
    pathway_id = paths$pathway_id,
    pathway_name = paths$pathway_name,
    global_pathway_size = m,
    CoreSeed_pathway_n = aobs,
    T_obs_logOR = Tobs,
    null_exceed_n = as.integer(exceed),
    B = Bn,
    P_emp = pemp,
    null_mean_T = null_mean,
    null_sd_T = sqrt(null_var)
  )
}

mdv4_add_qemp <- function(dt) {
  dt[, q_emp := p.adjust(P_emp, method = "BH"), by = source]
  dt[]
}

mdv4_supportive_direction <- function(s) {
  is.finite(s$Fisher_OR) & s$Fisher_OR > 1 &
    is.finite(s$MH_OR) & s$MH_OR > 1 &
    is.finite(s$goseq_over_P) & is.finite(s$goseq_under_P) &
    s$goseq_over_P < s$goseq_under_P
}
