# Per-study model fitting. Engines return the effect table (see dev/SPEC.md).

.mm_engines <- c("bezi_gamlss", "lm_asin", "lm_clr", "zibb_glmmtmb")
.mm_engine_scale <- c(bezi_gamlss = "logit_mean", lm_asin = "asin_sqrt_mean", lm_clr = "clr_mean",
                      zibb_glmmtmb = "logit_mean")
.zi_se_limit <- 10   # a logit-scale standard error above this means the effect is not identifiable

#' Fit a model to every feature of a study
#'
#' Fits the same model to each feature and returns the estimates as an effect
#' table. Four engines are available.
#'
#' * `"zibb_glmmtmb"`: a zero-inflated beta-binomial model on the counts with the
#'   sequencing depth as the number of trials (`glmmTMB`), a logit link for the mean
#'   abundance and an optional subject random intercept. Needs counts. Returns the
#'   abundance effects (`part = "abundance"`, scale `"logit_mean"`) and, if
#'   `prevalence` is given, effects on the probability of a structural zero
#'   (`part = "prevalence"`, scale `"logit_zero_prob"`; a positive value means the
#'   feature is more often absent).
#' * `"bezi_gamlss"`: the model of metamicrobiomeR 1.x: a GAMLSS with a zero-inflated
#'   beta (BEZI) family on relative abundance, logit link for the mean. Only the mean
#'   part depends on the covariates, as in 1.x. Scale `"logit_mean"`.
#' * `"lm_asin"`: a linear model on the arcsine square root of the relative abundance
#'   (the transformation used by MaAsLin), a linear mixed model with a subject random
#'   intercept when `longitudinal` (p-values from the normal approximation, as in 1.x),
#'   otherwise `lm()` with t-based p-values. Scale `"asin_sqrt_mean"`.
#' * `"lm_clr"`: the same linear (mixed) model on the centred log-ratio transform of the
#'   features at the analysed rank, as a compositional sensitivity analysis. A pseudocount
#'   replaces zeros. Scale `"clr_mean"`.
#'
#' Estimates from different engines are on different scales and are not comparable;
#' [bind_effects()] refuses to combine them.
#'
#' Features are kept if they are present (abundance above zero) in at least
#' `min_prevalence` of all samples and their mean relative abundance exceeds
#' `min_mean`. Samples with a missing value in the model are dropped for that feature.
#'
#' @param study An `mm_study` object.
#' @param formula One-sided formula with the covariates, e.g. `~ group + age`.
#'   Variables are looked up in the sample data.
#' @param engine Model engine (see Details), default `"bezi_gamlss"`.
#' @param rank Optional character vector of feature ranks to analyse (matched
#'   against the `rank` column of the feature data). `NULL` analyses all features.
#' @param longitudinal Add a subject-level random intercept? Defaults to `TRUE`
#'   when `study$subject` is set. The subject is the sample-data column named in
#'   `study$subject`.
#' @param min_prevalence,min_mean Feature filters (see Details).
#' @param n_cycles (`bezi_gamlss` only) Maximum number of outer iterations of the GAMLSS
#'   fitting algorithm (`n.cyc` in `gamlss::gamlss.control()`). The default is 100.
#'   Version 1.x used the GAMLSS default of 20, at which most fits on real data did not
#'   converge (see `dev/VALIDATION.md`); use `n_cycles = 20` to reproduce 1.x results.
#'   Fits that reach the limit are flagged `converged = FALSE` and a warning reports how many.
#' @param exclude_unassigned Drop features that are unassigned at their own rank
#'   (for example a genus-level feature with no genus, such as the 1.x
#'   `...f__lachnospiraceae.other`)? Such features lump together unrelated
#'   organisms. Default `TRUE`; use `FALSE` to reproduce 1.x.
#' @param prevalence (`zibb_glmmtmb` only) One-sided formula for the covariates of the
#'   probability of a structural zero, e.g. `~ group`. `NULL` (default) fits an
#'   intercept-only zero-inflation part, so only abundance effects are reported.
#'   Prevalence effects can be poorly identified; unidentifiable ones are flagged
#'   `converged = FALSE`.
#' @param zero_inflation (`zibb_glmmtmb` only) Include a zero-inflation part? It is
#'   dropped automatically for features without zero counts. `FALSE` fits a plain
#'   beta-binomial model.
#' @param pseudocount (`lm_clr` only) Value added to counts (default 0.5) or, for
#'   proportions, to relative abundances (default half the smallest positive value).
#' @return An `mm_effects` object with one row per feature, model part and covariate
#'   term. The attribute `fit_info` is a data.frame with one row per fitted feature
#'   (samples used, convergence, number of warnings, and any error message).
#' @export
fit_study <- function(study, formula, engine = "bezi_gamlss", rank = NULL,
                      longitudinal = NULL, min_prevalence = 0.05, min_mean = 5e-05,
                      n_cycles = 100, exclude_unassigned = TRUE, prevalence = NULL,
                      zero_inflation = TRUE, pseudocount = NULL) {
  if (!inherits(study, "mm_study")) stop("'study' must be an mm_study object.", call. = FALSE)
  engine <- match.arg(engine, .mm_engines)
  if (!inherits(formula, "formula") || length(formula) != 2) {
    stop("'formula' must be one-sided, e.g. ~ group + age.", call. = FALSE)
  }
  if (!is.null(prevalence) && (!inherits(prevalence, "formula") || length(prevalence) != 2)) {
    stop("'prevalence' must be a one-sided formula, e.g. ~ group.", call. = FALSE)
  }
  info <- mm_study_info(study)
  sdat <- mm_samples(study)
  vars <- unique(c(all.vars(formula), if (!is.null(prevalence)) all.vars(prevalence)))
  miss <- setdiff(vars, colnames(sdat))
  if (length(miss)) stop("Variable(s) not in the sample data: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  subject <- info$subject
  if (is.null(longitudinal)) longitudinal <- !is.na(subject)
  if (longitudinal && is.na(subject)) {
    stop("A longitudinal model needs study$subject (the repeated-measures ID column).", call. = FALSE)
  }
  counts <- mm_counts(study)
  if (engine == "zibb_glmmtmb" && is.null(counts)) {
    stop("The 'zibb_glmmtmb' engine needs counts; this study only has relative abundances.", call. = FALSE)
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
  composition <- idx                                   # the features that make up the composition (for CLR)
  n_all <- ncol(ra)
  keep <- vapply(idx, function(i) {
    x <- ra[i, ]
    sum(!is.na(x) & x > 0) >= min_prevalence * n_all && mean(x, na.rm = TRUE) > min_mean
  }, logical(1))
  idx <- idx[keep]

  Y <- switch(engine,
    bezi_gamlss = ra,
    zibb_glmmtmb = ra,
    lm_asin = asin(sqrt(ra)),
    lm_clr = .clr_matrix(study, composition, pseudocount))
  depth <- mm_depth(study)

  re_term <- if (engine == "bezi_gamlss") "random(.subject)" else "(1 | .subject)"
  terms_f <- attr(stats::terms(formula), "term.labels")
  f <- stats::reformulate(c(terms_f, if (longitudinal) re_term), response = ".y")

  rows <- vector("list", length(idx)); log <- vector("list", length(idx))
  for (k in seq_along(idx)) {
    i <- idx[k]; fid <- rownames(ra)[i]
    d <- d0
    d$.y <- Y[fid, ]
    if (engine == "zibb_glmmtmb") {
      d$.y <- counts[fid, ]; d$.fail <- depth - d$.y
      d <- d[depth > 0, , drop = FALSE]
    }
    if (engine == "bezi_gamlss") d$.y[which(d$.y >= 1)] <- 0.9999
    d <- stats::na.omit(d)
    res <- switch(engine,
      bezi_gamlss = .fit_bezi_gamlss(f, d, n_cycles),
      lm_asin = .fit_lm(f, d, longitudinal, .mm_engine_scale[["lm_asin"]]),
      lm_clr = .fit_lm(f, d, longitudinal, .mm_engine_scale[["lm_clr"]]),
      zibb_glmmtmb = .fit_zibb(terms_f, d, longitudinal, prevalence, zero_inflation))
    log[[k]] <- data.frame(feature = fid, n = nrow(d), converged = res$converged,
                           n_warnings = length(res$warnings), error = res$error, stringsAsFactors = FALSE)
    if (!is.null(res$rows) && nrow(res$rows)) {
      rk <- if ("rank" %in% colnames(feats)) feats$rank[i] else NA_character_
      rows[[k]] <- cbind(feature = fid, rank = rk, res$rows, n = nrow(d), stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, rows)
  if (is.null(tab)) tab <- data.frame(feature = character(), rank = character(), term = character(),
                                      part = character(), estimate = numeric(), se = numeric(), p_value = numeric(),
                                      df = numeric(), effect_scale = character(), converged = logical(), n = integer())
  out <- new_effect_table(study = rep(info$id, nrow(tab)), feature = tab$feature, term = tab$term,
                          estimate = tab$estimate, se = tab$se, rank = tab$rank, part = tab$part,
                          df = tab$df, p_value = tab$p_value, n = tab$n, engine = rep(engine, nrow(tab)),
                          effect_scale = tab$effect_scale, converged = tab$converged)
  fit_info <- do.call(rbind, log)
  attr(out, "fit_info") <- fit_info
  if (!is.null(fit_info) && any(!is.na(fit_info$error))) {
    warning(sum(!is.na(fit_info$error)), " of ", nrow(fit_info), " fits failed (first error: ",
            fit_info$error[!is.na(fit_info$error)][1], ").", call. = FALSE)
  }
  if (!is.null(fit_info) && any(!fit_info$converged)) {
    warning(sum(!fit_info$converged), " of ", nrow(fit_info), " fits did not converge",
            if (engine == "bezi_gamlss") paste0(" within ", n_cycles, " iterations"),
            " (see attr(, \"fit_info\")).", call. = FALSE)
  }
  out
}

# Centred log-ratio matrix (features x samples) over the features that make up the composition
.clr_matrix <- function(study, composition, pseudocount) {
  counts <- mm_counts(study); ra <- mm_relabund(study)
  X <- if (!is.null(counts)) counts[composition, , drop = FALSE] + (if (is.null(pseudocount)) 0.5 else pseudocount) else {
    r <- ra[composition, , drop = FALSE]
    r + (if (is.null(pseudocount)) min(r[r > 0], na.rm = TRUE) / 2 else pseudocount)
  }
  L <- log(X)
  out <- sweep(L, 2, colMeans(L), "-")
  out[, mm_depth(study) %in% 0 | colSums(is.na(ra)) > 0] <- NA_real_
  out
}

.coef_rows <- function(ct, part, scale, converged, p_col = "Pr(>|t|)") {
  ct <- ct[rownames(ct) != "(Intercept)", , drop = FALSE]
  if (!nrow(ct)) return(NULL)
  data.frame(term = rownames(ct), part = part, estimate = unname(ct[, "Estimate"]), se = unname(ct[, "Std. Error"]),
             p_value = unname(ct[, p_col]), df = NA_real_, effect_scale = scale, converged = converged,
             stringsAsFactors = FALSE)
}

# One GAMLSS-BEZI fit; mean-part coefficients without the intercept
.fit_bezi_gamlss <- function(f, d, n_cycles = 100) {
  warns <- character()
  catch <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  fit <- tryCatch(catch(gamlss::gamlss(f, family = "BEZI", data = d, trace = FALSE,
                                       control = gamlss::gamlss.control(n.cyc = n_cycles, trace = FALSE))),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(rows = NULL, converged = FALSE, warnings = warns, error = conditionMessage(fit)))
  }
  sm <- NULL
  utils::capture.output(sm <- tryCatch(catch(summary(fit, save = TRUE)), error = function(e) e))
  if (inherits(sm, "error")) {
    return(list(rows = NULL, converged = isTRUE(fit$converged), warnings = warns, error = conditionMessage(sm)))
  }
  k <- length(stats::coef(fit, "mu"))
  conv <- isTRUE(fit$converged)
  list(rows = .coef_rows(sm$coef.table[seq_len(k), , drop = FALSE], "abundance", .mm_engine_scale[["bezi_gamlss"]], conv),
       converged = conv, warnings = warns, error = NA_character_)
}

# One linear (mixed) model on a transformed scale
.fit_lm <- function(f, d, longitudinal, scale) {
  warns <- character()
  catch <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  fit <- tryCatch(
    if (longitudinal) catch(lme4::lmer(f, data = d)) else catch(stats::lm(f, data = d)),
    error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(rows = NULL, converged = FALSE, warnings = warns, error = conditionMessage(fit)))
  }
  ct <- catch(summary(fit))$coefficients
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
  rows <- .coef_rows(ct, "abundance", scale, conv)
  if (!is.null(rows)) rows$df <- df
  list(rows = rows, converged = conv, warnings = warns, error = NA_character_)
}

# One zero-inflated beta-binomial fit (glmmTMB)
.fit_zibb <- function(terms_f, d, longitudinal, prevalence, zero_inflation) {
  warns <- character()
  catch <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  use_zi <- isTRUE(zero_inflation) && any(d$.y == 0)
  f <- stats::reformulate(c(terms_f, if (longitudinal) "(1 | .subject)"), response = "cbind(.y, .fail)")
  zif <- if (!use_zi) ~0 else if (is.null(prevalence)) ~1 else prevalence
  fit <- tryCatch(catch(glmmTMB::glmmTMB(f, ziformula = zif, family = glmmTMB::betabinomial(), data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(rows = NULL, converged = FALSE, warnings = warns, error = conditionMessage(fit)))
  }
  conv <- isTRUE(fit$fit$convergence == 0) && isTRUE(fit$sdr$pdHess)
  sm <- tryCatch(catch(summary(fit))$coefficients, error = function(e) e)
  if (inherits(sm, "error")) {
    return(list(rows = NULL, converged = conv, warnings = warns, error = conditionMessage(sm)))
  }
  rows <- .coef_rows(sm$cond, "abundance", .mm_engine_scale[["zibb_glmmtmb"]], conv, "Pr(>|z|)")
  if (use_zi && !is.null(prevalence) && !is.null(sm$zi)) {
    zr <- .coef_rows(sm$zi, "prevalence", "logit_zero_prob", conv, "Pr(>|z|)")
    if (!is.null(zr)) {
      zr$converged <- conv & is.finite(zr$se) & zr$se <= .zi_se_limit
      rows <- rbind(rows, zr)
    }
  }
  list(rows = rows, converged = conv, warnings = warns, error = NA_character_)
}
