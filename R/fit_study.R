# Per-study model fitting. Engines return the effect table (see dev/SPEC.md).

.mm_engines <- c("bezi_gamlss")

#' Fit a model to every feature of a study
#'
#' Fits the same model to each feature and returns the estimates as an effect
#' table. The first engine, `"bezi_gamlss"`, is the model of metamicrobiomeR 1.x:
#' a GAMLSS with a zero-inflated beta (BEZI) family on relative abundance, with a
#' logit link for the mean. Only the mean (abundance) part depends on the
#' covariates; the zero-inflation and dispersion parts are intercept-only, as in 1.x.
#'
#' Features are kept if they are present (abundance above zero) in at least
#' `min_prevalence` of all samples and their mean relative abundance exceeds
#' `min_mean`. Relative abundances of 1 or more are set to 0.9999 and samples with
#' missing values in the model are dropped for that feature.
#'
#' @param study An `mm_study` object.
#' @param formula One-sided formula with the covariates, e.g. `~ group + age`.
#'   Variables are looked up in the sample data.
#' @param engine Model engine. Currently only `"bezi_gamlss"`.
#' @param rank Optional character vector of feature ranks to analyse (matched
#'   against the `rank` column of the feature data). `NULL` analyses all features.
#' @param longitudinal Add a subject-level random intercept? Defaults to `TRUE`
#'   when `study$subject` is set. The subject is the sample-data column named in
#'   `study$subject`.
#' @param min_prevalence,min_mean Feature filters (see Details).
#' @param n_cycles Maximum number of outer iterations of the GAMLSS fitting
#'   algorithm (`n.cyc` in `gamlss::gamlss.control()`). The default is 100.
#'   Version 1.x used the GAMLSS default of 20, at which most fits on real data
#'   did not converge (see `dev/VALIDATION.md`); use `n_cycles = 20` to reproduce
#'   1.x results. Fits that reach the limit are flagged `converged = FALSE` and a
#'   warning reports how many.
#' @return An `mm_effects` object with one row per feature and covariate term,
#'   `part = "abundance"` and `effect_scale = "logit_mean"`. The attribute
#'   `fit_info` is a data.frame with one row per fitted feature (samples used,
#'   convergence, number of warnings, and any error message).
#' @export
fit_study <- function(study, formula, engine = "bezi_gamlss", rank = NULL,
                      longitudinal = NULL, min_prevalence = 0.05, min_mean = 5e-05,
                      n_cycles = 100) {
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
  n_all <- ncol(ra)
  keep <- vapply(idx, function(i) {
    x <- ra[i, ]
    sum(!is.na(x) & x > 0) >= min_prevalence * n_all && mean(x, na.rm = TRUE) > min_mean
  }, logical(1))
  idx <- idx[keep]

  rhs <- c(attr(stats::terms(formula), "term.labels"), if (longitudinal) "random(.subject)")
  f <- stats::reformulate(rhs, response = ".y")

  rows <- vector("list", length(idx)); log <- vector("list", length(idx))
  for (k in seq_along(idx)) {
    i <- idx[k]
    d <- d0
    d$.y <- ra[i, ]
    d$.y[which(d$.y >= 1)] <- 0.9999
    d <- stats::na.omit(d)
    res <- .fit_bezi_gamlss(f, d, n_cycles)
    fid <- rownames(ra)[i]
    log[[k]] <- data.frame(feature = fid, n = nrow(d), converged = res$converged,
                           n_warnings = length(res$warnings), error = res$error,
                           stringsAsFactors = FALSE)
    if (!is.null(res$table) && nrow(res$table)) {
      rk <- if ("rank" %in% colnames(feats)) feats$rank[i] else NA_character_
      rows[[k]] <- data.frame(feature = fid, rank = rk, term = rownames(res$table),
                              estimate = res$table[, "Estimate"], se = res$table[, "Std. Error"],
                              p_value = res$table[, "Pr(>|t|)"], n = nrow(d),
                              converged = res$converged, stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, rows)
  if (is.null(tab)) tab <- data.frame(feature = character(), rank = character(), term = character(),
                                      estimate = numeric(), se = numeric(), p_value = numeric(),
                                      n = integer(), converged = logical())
  out <- new_effect_table(study = info$id, feature = tab$feature, term = tab$term,
                          estimate = tab$estimate, se = tab$se, rank = tab$rank,
                          part = "abundance", p_value = tab$p_value, n = tab$n,
                          engine = engine, effect_scale = "logit_mean", converged = tab$converged)
  fit_info <- do.call(rbind, log)
  attr(out, "fit_info") <- fit_info
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
  list(table = ct, converged = isTRUE(fit$converged), warnings = warns, error = NA_character_)
}
