# Shared numerical-certification solver for MDV-3 runtime patch v1p3.
# Scientific model/estimator target is unchanged. AS_mean remains first-line.
# MPL_Jeffreys(a=1/2) is a numerical-equivalent fallback for logistic regression.

mdv3_check_brglm2_version <- function(cfg) {
  minv <- as.character(cfg$numerical_certification$require_brglm2_min_version)
  cur <- as.character(utils::packageVersion("brglm2"))
  if (utils::compareVersion(cur, minv) < 0) {
    stop("brglm2 >= ", minv, " is required for the numerical-certification patch; installed=", cur)
  }
  invisible(cur)
}

mdv3_safe_exp <- function(x) {
  hi <- log(.Machine$double.xmax)
  lo <- log(.Machine$double.xmin)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  out[ok & x > hi] <- Inf
  out[ok & x < lo] <- 0
  mid <- ok & x >= lo & x <= hi
  out[mid] <- exp(x[mid])
  out
}

mdv3_warning_flags <- function(warnings) {
  txt <- paste(unique(warnings), collapse = " | ")
  list(
    text = txt,
    nonconvergence = grepl("algorithm did not converge", txt, fixed = TRUE),
    prob01 = grepl("fitted probabilities numerically 0 or 1 occurred", txt, fixed = TRUE),
    boundary = grepl("algorithm stopped at boundary value", txt, fixed = TRUE)
  )
}

mdv3_attempt_control <- function(a) {
  brglm2::brglm_control(
    epsilon = as.numeric(a$epsilon),
    maxit = as.integer(a$maxit),
    type = as.character(a$type),
    slowit = as.numeric(a$slowit),
    response_adjustment = as.numeric(a$response_adjustment),
    max_step_factor = as.integer(a$max_step_factor),
    a = as.numeric(a$a),
    trace = FALSE,
    check_aliasing = TRUE
  )
}

mdv3_fit_once <- function(formula, data, attempt) {
  warnings <- character()
  ctrl <- mdv3_attempt_control(attempt)
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm(
        formula,
        data = data,
        family = stats::binomial("logit"),
        method = brglm2::brglmFit,
        control = ctrl
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(fit = fit, warnings = unique(warnings), flags = mdv3_warning_flags(warnings))
}

mdv3_certified_fit <- function(formula, data, cfg, coef_name = "pathway_member", ci_level = 0.95) {
  attempts <- cfg$numerical_certification$attempts
  zcrit <- stats::qnorm(1 - (1 - as.numeric(ci_level)) / 2)
  attempt_log <- list()

  for (i in seq_along(attempts)) {
    a <- attempts[[i]]
    one <- mdv3_fit_once(formula, data, a)
    fit <- one$fit
    rec <- list(
      attempt = as.character(a$id),
      solver_type = as.character(a$type),
      converged = FALSE,
      iterations = NA_integer_,
      warning = one$flags$text,
      beta = NA_real_,
      se = NA_real_
    )

    if (!inherits(fit, "error")) {
      rec$converged <- isTRUE(fit$converged)
      rec$iterations <- suppressWarnings(as.integer(fit$iter))
      cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
      vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
      if (!is.null(cf) && !is.null(vc) && coef_name %in% names(cf) && coef_name %in% rownames(vc)) {
        rec$beta <- as.numeric(cf[[coef_name]])
        rec$se <- sqrt(as.numeric(vc[coef_name, coef_name]))
      }
    } else {
      rec$warning <- paste(c(rec$warning, paste0("ERROR: ", conditionMessage(fit))), collapse = " | ")
    }
    attempt_log[[length(attempt_log) + 1L]] <- rec

    valid <- isTRUE(rec$converged) && is.finite(rec$beta) && is.finite(rec$se) && rec$se > 0 && !one$flags$nonconvergence && !one$flags$boundary
    if (valid) {
      beta <- rec$beta
      se <- rec$se
      p <- 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE)
      log_lo <- beta - zcrit * se
      log_hi <- beta + zcrit * se
      return(list(
        certified = TRUE,
        fit = fit,
        beta = beta,
        se = se,
        p = p,
        log_ci_low = log_lo,
        log_ci_high = log_hi,
        or = mdv3_safe_exp(beta),
        ci_low = mdv3_safe_exp(log_lo),
        ci_high = mdv3_safe_exp(log_hi),
        solver_attempt = rec$attempt,
        solver_type = rec$solver_type,
        converged = TRUE,
        iterations = rec$iterations,
        warning = one$flags$text,
        prob01_warning = one$flags$prob01,
        attempt_log = attempt_log
      ))
    }
  }

  list(
    certified = FALSE,
    fit = NULL,
    beta = NA_real_, se = NA_real_, p = 1,
    log_ci_low = NA_real_, log_ci_high = NA_real_,
    or = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    solver_attempt = paste(vapply(attempt_log, `[[`, character(1), "attempt"), collapse = ">"),
    solver_type = paste(vapply(attempt_log, `[[`, character(1), "solver_type"), collapse = ">"),
    converged = FALSE,
    iterations = suppressWarnings(max(vapply(attempt_log, function(x) as.integer(if (is.na(x$iterations)) 0L else x$iterations), integer(1)))),
    warning = paste(vapply(attempt_log, `[[`, character(1), "warning"), collapse = " || "),
    prob01_warning = any(vapply(attempt_log, function(x) grepl("fitted probabilities numerically 0 or 1 occurred", x$warning, fixed = TRUE), logical(1))),
    attempt_log = attempt_log
  )
}
