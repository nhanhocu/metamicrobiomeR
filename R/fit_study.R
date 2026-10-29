# Per-study model fitting. Engines return the effect table (see dev/SPEC.md).

.mm_engines <- c("bezi_gamlss", "lm_asin")

#' Fit a model to every feature of a study
#'
#' Fits the same model to each feature and returns the estimates as an effect
#' table. Two engines from metamicrobiomeR 1.x are available.
#'
#' * `"bezi_gamlss"`: a GAMLSS with a zero-inflated beta (BEZI) family on relative
#'   abundance, with a logit link for the mean. Only the mean (abundance) part
#'   depends on the covariates; the zero-inflation and dispersion parts are
#'   intercept-only, as in 1.x. Effects are on the logit scale (`"logit_mean"`).
#' * `"lm_asin"`: a linear model on the arcsine square root of the relative
#'   abundance (the transformation used by MaAsLin), a linear mixed model with a
#'   subject random intercept when `longitudinal` (p-values from the normal
#'   approximation, as in 1.x), otherwise `lm()` with t-based p-values. Effects are
#'   on the arcsine-square-root scale (`"asin_sqrt_mean"`) and are not comparable
#'   with `"bezi_gamlss"` estimates.
#'
#' Features are kept if they are present (abundance above zero) in at least
#' `min_prevalence` of all samples and their mean relative abundance exceeds
#' `min_mean`. Relative abundances of 1 or more are set to 0.9999 and samples with
#' missing values in the model are dropped for that feature.
#'
#' @param study An `mm_study` object.
#' @param formula One-sided formula with the covariates, e.g. `~ group + age`.
#'   Variables are looked up in the sample data.
#' @param engine Model engine: `"bezi_gamlss"` (default) or `"lm_asin"`.
#' @param rank Optional character vector of feature ranks to analyse (matched
#'   against the `rank` column of the feature data). `NULL` analyses all features.
#' @param longitudinal Add a subject-level random intercept? Defaults to `TRUE`
#'   when `study$subject` is set. The subject is the sample-data column named in
#'   `study$subject`.
#' @param min_prevalence,min_mean Feature filters (see Details).
#' @param exclude_unassigned Drop features that are unassigned at their own rank
#'   (for example a genus-level feature with no genus, such as the 1.x
#'   `...f__lachnospiraceae.other`)? Such features lump together unrelated
#'   organisms. Default `TRUE`; use `FALSE` to reproduce 1.x.
#' @param n_cycles (`bezi_gamlss` only) Maximum number of outer iterations of the GAMLSS fitting
#'   algorithm (`n.cyc` in `gamlss::gamlss.control()`). The default is 100.
#'   Version 1.x used the GAMLSS default of 20, at which most fits on real data
#'   did not converge (see `dev/VALIDATION.md`); use `n_cycles = 20` to reproduce
#'   1.x results. Fits that reach the limit are flagged `converged = FALSE` and a
#'   warning reports how many.
#' @return An `mm_effects` object with one row per feature and covariate term,
#'   `part = "abundance"` and the effect scale of the engine. The attribute
#'   `fit_info` is a data.frame with one row per fitted feature (samples used,
#'   convergence, number of warnings, and any error message).
#' @export
fit_study <- function(study, formula, engine = "bezi_gamlss", rank = NULL,
                      longitudinal = NULL, min_prevalence = 0.05, min_mean = 5e-05,
                      n_cycles = 100, exclude_unassigned = TRUE) {
  if (!inherits(study, "mm_study")) stop("'study' must be an mm_study object.", call. = FALSE)
  engine <- match.arg(engine, .mm_engines)
  if (!inherits(formula, "formula") || length(formula) != 2) {
    stop("'formula' must be one-sided, e.g. ~ group + age.", call. = FALSE)
  }
  info <- mm_study_info(study)
  sdat <- mm_samples(study)
  vars <- all.vars(formula)
  miss <- setdiff(vars, colnames(sdat))
  if (length(miss)) stop("Variable(s) not in the sample data: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  subject <- info$subject
  if (is.null(longitudinal)) longitudinal <- !is.na(subject)
  if (longitudinal && is.na(subject)) {
    stop("A longitudinal model needs study$subject (the repeated-measures ID column).", call. = FALSE)
  }

  d0 <- sdat[, vars, drop = FALSE]
  for (v in vars) if (is.factor(d0[[v]])) d0[[v]] <- droplevels(d0[[v]])
  if (longitudinal) d0$.subject <- factor(sdat[[subject]])

  ra <- mm_relabund(study)
  feats <- mm_features(study)
  idx <- seq_len(nrow(ra))
  if (!is.null(rank)) {
    if (!"rank" %in% colnames(feats)) stop("The feature data has no 'rank' column.", call. = FALSE)
    idx <- which(feats$rank %in% rank)
    if (!length(idx)) stop("No features with rank: ", paste(rank, collapse = ", "), ".", call. = FALSE)
  }
  if (exclude_unassigned && "rank" %in% colnames(feats)) {
    unassigned <- vapply(idx, function(i) {
      r <- feats$rank[i]
      !is.na(r) && r %in% colnames(feats) && is.na(feats[[r]][i])
    }, logical(1))
    idx <- idx[!unassigned]
  }
  n_all <- ncol(ra)
  keep <- vapply(idx, function(i) {
    x <- ra[i, ]
    sum(!is.na(x) & x > 0) >= min_prevalence * n_all && mean(x, na.rm = TRUE) > min_mean
  }, logical(1))
  idx <- idx[keep]

  re_term <- if (engine == "bezi_gamlss") "random(.subject)" else "(1 | .subject)"
  rhs <- c(attr(stats::terms(formula), "term.labels"), if (longitudinal) re_term)
  f <- stats::reformulate(rhs, response = ".y")

  rows <- vector("list", length(idx)); log <- vector("list", length(idx))
  for (k in seq_along(idx)) {
    i <- idx[k]
    d <- d0
    d$.y <- ra[i, ]
    if (engine == "bezi_gamlss") d$.y[which(d$.y >= 1)] <- 0.9999 else d$.y <- asin(sqrt(d$.y))
    d <- stats::na.omit(d)
    res <- switch(engine,
                  bezi_gamlss = .fit_bezi_gamlss(f, d, n_cycles),
                  lm_asin = .fit_lm_asin(f, d, longitudinal))
    fid <- rownames(ra)[i]
    log[[k]] <- data.frame(feature = fid, n = nrow(d), converged = res$converged,
                           n_warnings = length(res$warnings), error = res$error,
                           stringsAsFactors = FALSE)
    if (!is.null(res$table) && nrow(res$table)) {
      rk <- if ("rank" %in% colnames(feats)) feats$rank[i] else NA_character_
      rows[[k]] <- data.frame(feature = fid, rank = rk, term = rownames(res$table),
                              estimate = res$table[, "Estimate"], se = res$table[, "Std. Error"],
                              p_value = res$table[, "Pr(>|t|)"], df = res$df, n = nrow(d),
                              converged = res$converged, stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, rows)
  if (is.null(tab)) tab <- data.frame(feature = character(), rank = character(), term = character(),
                                      estimate = numeric(), se = numeric(), p_value = numeric(), df = numeric(),
                                      n = integer(), converged = logical())
  nr <- nrow(tab)
  out <- new_effect_table(study = rep(info$id, nr), feature = tab$feature, term = tab$term,
                          estimate = tab$estimate, se = tab$se, rank = tab$rank,
                          part = rep("abundance", nr), df = tab$df, p_value = tab$p_value, n = tab$n,
                          engine = rep(engine, nr),
                          effect_scale = rep(if (engine == "bezi_gamlss") "logit_mean" else "asin_sqrt_mean", nr),
                          converged = tab$converged)
  fit_info <- do.call(rbind, log)
  attr(out, "fit_info") <- fit_info
  if (!is.null(fit_info) && any(!is.na(fit_info$error))) {
    warning(sum(!is.na(fit_info$error)), " of ", nrow(fit_info), " fits failed (first error: ",
            fit_info$error[!is.na(fit_info$error)][1], ").", call. = FALSE)
  }
  if (!is.null(fit_info) && any(!fit_info$converged)) {
    warning(sum(!fit_info$converged), " of ", nrow(fit_info), " fits did not converge within ", n_cycles,
            " iterations (see attr(, \"fit_info\")); consider a larger 'n_cycles'.", call. = FALSE)
  }
  out
}

# One GAMLSS-BEZI fit; returns the mean-part coefficient table without the intercept
.fit_bezi_gamlss <- function(f, d, n_cycles = 100) {
  warns <- character()
  catch <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  fit <- tryCatch(catch(gamlss::gamlss(f, family = "BEZI", data = d, trace = FALSE,
                                                 control = gamlss::gamlss.control(n.cyc = n_cycles, trace = FALSE))),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(table = NULL, converged = FALSE, warnings = warns, error = conditionMessage(fit)))
  }
  sm <- NULL
  utils::capture.output(sm <- tryCatch(catch(summary(fit, save = TRUE)), error = function(e) e))
  if (inherits(sm, "error")) {
    return(list(table = NULL, converged = isTRUE(fit$converged), warnings = warns, error = conditionMessage(sm)))
  }
  k <- length(stats::coef(fit, "mu"))
  ct <- sm$coef.table[seq_len(k), , drop = FALSE]
  ct <- ct[rownames(ct) != "(Intercept)", , drop = FALSE]
  list(table = ct, converged = isTRUE(fit$converged), warnings = warns, error = NA_character_, df = NA_real_)
}

# One linear (mixed) model on the arcsine square root scale
.fit_lm_asin <- function(f, d, longitudinal) {
  warns <- character()
  catch <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  fit <- tryCatch(
    if (longitudinal) catch(lme4::lmer(f, data = d)) else catch(stats::lm(f, data = d)),
    error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(table = NULL, converged = FALSE, warnings = warns, error = conditionMessage(fit), df = NA_real_))
  }
  sm <- catch(summary(fit))
  ct <- sm$coefficients
  if (longitudinal) {
    # lme4 gives no p-values: normal approximation, as in 1.x
    ct <- cbind(ct[, c("Estimate", "Std. Error", "t value"), drop = FALSE],
                `Pr(>|t|)` = 2 * stats::pnorm(-abs(ct[, "t value"])))
    conv <- isTRUE(fit@optinfo$conv$opt == 0) && !any(grepl("converge", warns))
    df <- NA_real_
  } else {
    conv <- TRUE
    df <- as.numeric(stats::df.residual(fit))
  }
  ct <- ct[rownames(ct) != "(Intercept)", , drop = FALSE]
  list(table = ct, converged = conv, warnings = warns, error = NA_character_, df = df)
}
