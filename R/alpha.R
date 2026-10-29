# Alpha diversity: per-sample diversity indices and their group effects

#' Per-sample alpha diversity
#'
#' @param study An `mm_study`.
#' @param metrics Any of `"shannon"` (natural log), `"simpson"` (Gini-Simpson, 1 - sum p^2),
#'   `"inv_simpson"`, `"pielou"` (Shannon evenness), `"observed"` (features with a count above
#'   zero) and `"chao1"` (bias-corrected). The last two need counts.
#' @param rarefy Optional integer depth to subsample every sample to (without replacement)
#'   before computing the indices; samples with fewer reads are dropped. Richness-type indices
#'   depend strongly on sequencing depth, so rarefying or a depth covariate is advisable.
#' @param seed Optional seed for the rarefaction.
#' @return A data.frame with one row per sample (row names are the sample IDs) and one column per metric.
#' @export
alpha_diversity <- function(study, metrics = c("shannon", "simpson", "observed"), rarefy = NULL, seed = NULL) {
  if (!inherits(study, "mm_study")) stop("'study' must be an mm_study object.", call. = FALSE)
  metrics <- match.arg(metrics, c("shannon", "simpson", "inv_simpson", "pielou", "observed", "chao1"), several.ok = TRUE)
  cnt <- mm_counts(study)
  if (any(metrics %in% c("observed", "chao1")) && is.null(cnt)) stop("'observed' and 'chao1' need counts.", call. = FALSE)
  if (!is.null(rarefy)) {
    if (is.null(cnt)) stop("Rarefaction needs counts.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)
    ok <- colSums(cnt) >= rarefy
    if (!any(ok)) stop("No sample has at least ", rarefy, " reads.", call. = FALSE)
    cnt <- cnt[, ok, drop = FALSE]
    cnt <- apply(cnt, 2, function(x) tabulate(sample(rep.int(seq_along(x), x), rarefy), nbins = length(x)))
    rownames(cnt) <- rownames(mm_counts(study))
  }
  p <- if (is.null(cnt)) mm_relabund(study) else sweep(cnt, 2, colSums(cnt), "/")
  samples <- colnames(p)
  H <- -colSums(ifelse(p > 0, p * log(p), 0))
  out <- data.frame(row.names = samples)
  for (m in metrics) out[[m]] <- switch(m,
    shannon = H, simpson = 1 - colSums(p^2), inv_simpson = 1 / colSums(p^2),
    observed = colSums(cnt > 0),
    pielou = ifelse(colSums(p > 0) > 1, H / log(colSums(p > 0)), NA_real_),
    chao1 = apply(cnt, 2, function(x) { f1 <- sum(x == 1); f2 <- sum(x == 2); sum(x > 0) + f1 * (f1 - 1) / (2 * (f2 + 1)) }))
  out
}

#' Group effects on alpha diversity
#'
#' Fits a linear model (or a linear mixed model with a subject random intercept for repeated
#' measures) to each alpha-diversity metric and returns the effects as an effect table, so
#' they can be pooled across studies with [pool_effects()] like any other feature.
#' By default each metric is standardised to mean 0 and standard deviation 1 within the
#' study first, which makes effects comparable between studies that measure on different
#' scales (as in 1.x); the estimates are then in standard-deviation units.
#'
#' @param study An `mm_study`.
#' @param formula One-sided formula with the covariates, e.g. `~ group + age`.
#' @param metrics,rarefy,seed As in [alpha_diversity()].
#' @param standardise Standardise each metric within the study? Default `TRUE`.
#' @param longitudinal Add a subject random intercept? Defaults to `TRUE` when `study$subject` is set.
#' @return An `mm_effects` object with `feature` = metric, `part = "overall"` and
#'   `effect_scale` `"sd_units"` (standardised) or `"alpha_raw"`.
#' @export
fit_alpha <- function(study, formula, metrics = c("shannon", "simpson", "observed"), standardise = TRUE,
                      longitudinal = NULL, rarefy = NULL, seed = NULL) {
  if (!inherits(study, "mm_study")) stop("'study' must be an mm_study object.", call. = FALSE)
  if (!inherits(formula, "formula") || length(formula) != 2) stop("'formula' must be one-sided, e.g. ~ group + age.", call. = FALSE)
  info <- mm_study_info(study); sdat <- mm_samples(study)
  vars <- all.vars(formula); miss <- setdiff(vars, colnames(sdat))
  if (length(miss)) stop("Variable(s) not in the sample data: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  if (is.null(longitudinal)) longitudinal <- !is.na(info$subject)
  if (longitudinal && is.na(info$subject)) stop("A longitudinal model needs study$subject.", call. = FALSE)
  al <- alpha_diversity(study, metrics, rarefy, seed)
  d0 <- sdat[rownames(al), vars, drop = FALSE]
  for (v in vars) if (is.factor(d0[[v]])) d0[[v]] <- droplevels(d0[[v]])
  if (longitudinal) d0$.subject <- factor(sdat[rownames(al), info$subject])
  f <- stats::reformulate(c(attr(stats::terms(formula), "term.labels"), if (longitudinal) "(1 | .subject)"), response = ".y")
  scale <- if (standardise) "sd_units" else "alpha_raw"
  rows <- list(); log <- list()
  for (m in colnames(al)) {
    d <- d0; y <- al[[m]]
    if (standardise) y <- (y - mean(y, na.rm = TRUE)) / stats::sd(y, na.rm = TRUE)
    d$.y <- y; d <- stats::na.omit(d)
    res <- .fit_lm(f, d, longitudinal, scale)
    log[[m]] <- data.frame(feature = m, n = nrow(d), converged = res$converged, n_warnings = length(res$warnings), error = res$error)
    if (!is.null(res$rows)) rows[[m]] <- cbind(feature = m, rank = NA_character_, res$rows, n = nrow(d), stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows)
  if (is.null(tab)) stop("No metric could be fitted.", call. = FALSE)
  tab$part <- "overall"
  out <- new_effect_table(study = rep(info$id, nrow(tab)), feature = tab$feature, term = tab$term, estimate = tab$estimate,
                          se = tab$se, rank = tab$rank, part = tab$part, df = tab$df, p_value = tab$p_value, n = tab$n,
                          engine = rep("lm_alpha", nrow(tab)), effect_scale = tab$effect_scale, converged = tab$converged)
  attr(out, "fit_info") <- do.call(rbind, log)
  out
}
