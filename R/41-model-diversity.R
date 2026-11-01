# Alpha diversity and beta diversity (PERMANOVA) per study, and combining beta diversity

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
#' @family per-study models
#' @export
alpha_diversity <- function(study, metrics = c("shannon", "simpson", "observed"), rarefy = NULL, seed = NULL) {
  .check_class(study, "mm_study", "study", "as_study()")
  metrics <- match.arg(metrics, c("shannon", "simpson", "inv_simpson", "pielou", "observed", "chao1"), several.ok = TRUE)
  cnt <- mm_counts(study)
  if (any(metrics %in% c("observed", "chao1")) && is.null(cnt)) stop("'observed' and 'chao1' need counts.", call. = FALSE)
  if (!is.null(rarefy)) {
    if (is.null(cnt)) stop("Rarefaction needs counts.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)
    ok <- colSums(cnt) >= rarefy
    if (!any(ok)) stop("No sample has at least ", format(rarefy, big.mark = ","), " reads; the deepest has ", format(max(colSums(cnt)), big.mark = ","), ".", call. = FALSE)
    if (!all(ok)) .inform("Rarefying to ", format(rarefy, big.mark = ","), " reads: ", sum(!ok), " of ", length(ok), " samples with fewer reads are dropped.")
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
#' scales; the estimates are then in standard-deviation units.
#'
#' @param study An `mm_study`.
#' @param formula One-sided formula with the covariates, e.g. `~ group + age`.
#' @param metrics,rarefy,seed As in [alpha_diversity()].
#' @param standardise Standardise each metric within the study? Default `TRUE`.
#' @param longitudinal Add a subject random intercept? Defaults to `TRUE` when `study$subject` is set.
#' @return An `mm_effects` object with `feature` = metric, `part = "overall"` and
#'   `effect_scale` `"sd_units"` (standardised) or `"alpha_raw"`.
#' @family per-study models
#' @export
fit_alpha <- function(study, formula, metrics = c("shannon", "simpson", "observed"), standardise = TRUE,
                      longitudinal = NULL, rarefy = NULL, seed = NULL) {
  .check_class(study, "mm_study", "study", "as_study()")
  if (!inherits(formula, "formula") || length(formula) != 2) stop("'formula' must be one-sided, e.g. ~ group + age.", call. = FALSE)
  info <- mm_study_info(study); sdat <- mm_samples(study)
  vars <- all.vars(formula); miss <- setdiff(vars, colnames(sdat))
  if (length(miss)) .missing_vars(miss, colnames(sdat))
  if (is.null(longitudinal)) longitudinal <- !is.na(info$subject)
  if (longitudinal && is.na(info$subject)) stop("A longitudinal model needs study$subject, the sample column that identifies the subject (see as_study()).", call. = FALSE)
  if (longitudinal) .need("lme4", "a subject random effect (longitudinal = FALSE avoids it)")
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

#' Per-study PERMANOVA
#'
#' Tests whether the covariates explain differences in community composition between
#' samples, with `vegan::adonis2()` on a distance matrix of the relative abundances. Each term
#' is tested marginally (given the others). For repeated measures (a subject is set) whole
#' subjects are permuted, which is valid only for terms that are constant within a subject; use
#' it for between-subject factors.
#'
#' @param study An `mm_study`.
#' @param formula One-sided formula with the covariates, e.g. `~ group`.
#' @param distance A `vegan::vegdist()` method (e.g. `"bray"`, `"jaccard"`) or `"aitchison"`
#'   (Euclidean distance on the centred log-ratio with a pseudocount).
#' @param permutations Number of permutations.
#' @param rank Optional ranks to include (as in [fit_study()]); `NULL` uses all features.
#' @param seed Optional seed for the permutations.
#' @param longitudinal Permute whole subjects? Defaults to `TRUE` when `study$subject` is set.
#'   Then every term must be constant within subjects, otherwise the function stops: such a term
#'   cannot be tested by permuting subjects. `FALSE` ignores the repeated measures, which makes
#'   the p-values too small.
#' @return A data.frame of class `mm_beta` with one row per term: `study`, `term`, `df`, `R2`
#'   (the share of the distance variation explained), `F`, `p_value`, `n` and `distance`.
#' @family per-study models
#' @export
fit_beta <- function(study, formula, distance = "bray", permutations = 999, rank = NULL, seed = NULL,
                     longitudinal = NULL) {
  .need("vegan", "fit_beta()"); .need("permute", "fit_beta()")
  .check_class(study, "mm_study", "study", "as_study()")
  if (!inherits(formula, "formula") || length(formula) != 2) stop("'formula' must be one-sided, e.g. ~ group.", call. = FALSE)
  info <- mm_study_info(study); sdat <- mm_samples(study)
  vars <- all.vars(formula); miss <- setdiff(vars, colnames(sdat))
  if (length(miss)) .missing_vars(miss, colnames(sdat))
  ra <- mm_relabund(study); feats <- mm_features(study)
  if (!is.null(rank)) {
    if (!"rank" %in% colnames(feats)) stop("The feature data has no 'rank' column, so 'rank' cannot select features; agglomerate() creates one.", call. = FALSE)
    ra <- ra[feats$rank %in% rank, , drop = FALSE]
    if (!nrow(ra)) stop("No features with rank ", paste0("'", rank, "'", collapse = ", "), "; the study has: ", .show(sort(unique(feats$rank)), 6), ".", call. = FALSE)
  }
  d0 <- sdat[, vars, drop = FALSE]
  ok <- stats::complete.cases(d0) & colSums(is.na(ra)) == 0 & colSums(ra) > 0
  if (sum(ok) < 4) stop("Too few usable samples: ", sum(ok), " of ", length(ok), " have complete model variables and some abundance, but at least 4 are needed.", call. = FALSE)
  if (!all(ok)) .inform(sum(!ok), " of ", length(ok), " samples are left out (missing model variables or no abundance).")
  d0 <- d0[ok, , drop = FALSE]; ra <- ra[, ok, drop = FALSE]
  for (v in vars) if (is.factor(d0[[v]])) d0[[v]] <- droplevels(d0[[v]])
  X <- t(ra)
  dm <- if (distance == "aitchison") {
    pc <- min(X[X > 0]) / 2; L <- log(X + pc); stats::dist(L - rowMeans(L))
  } else vegan::vegdist(X, method = distance)
  if (is.null(longitudinal)) longitudinal <- !is.na(info$subject)
  if (longitudinal && is.na(info$subject)) stop("A longitudinal analysis needs study$subject.", call. = FALSE)
  subject <- if (longitudinal) factor(sdat[ok, info$subject]) else NULL
  if (longitudinal) {
    varying <- vars[vapply(vars, function(v) any(tapply(d0[[v]], subject, function(x) length(unique(x))) > 1), logical(1))]
    if (length(varying)) stop(paste0("'", varying, "'", collapse = ", "), if (length(varying) == 1) " varies" else " vary", " within subjects. ",
                              "Permuting whole subjects can only test between-subject terms. ",
                              "Use longitudinal = FALSE to ignore the repeated measures (p-values will be too small).", call. = FALSE)
  }
  perm <- if (is.null(subject)) permutations else
    permute::how(nperm = permutations, plots = permute::Plots(strata = subject, type = "free"), within = permute::Within(type = "none"))
  if (!is.null(seed)) set.seed(seed)
  rhs <- paste(attr(stats::terms(formula), "term.labels"), collapse = " + ")
  fit <- vegan::adonis2(stats::as.formula(paste("dm ~", rhs)), data = d0, permutations = perm, by = "margin")
  tab <- as.data.frame(fit); tab <- tab[!rownames(tab) %in% c("Residual", "Total"), , drop = FALSE]
  out <- data.frame(study = info$id, term = rownames(tab), df = tab$Df, R2 = tab$R2, F = tab$F, p_value = tab[["Pr(>F)"]],
                    n = nrow(d0), distance = distance, stringsAsFactors = FALSE)
  class(out) <- c("mm_beta", "data.frame")
  out
}

#' Combine PERMANOVA results across studies
#'
#' PERMANOVA gives no estimate with a standard error, so studies cannot be pooled like effect
#' sizes. This combines the permutation p-values of the studies for each term, either with
#' Fisher's method or with weighted Stouffer's method (weights proportional to the square root of
#' the sample size), and reports the sample-size-weighted mean R-squared as a descriptive
#' summary. Heterogeneity is not modelled, and the studies must use the same distance.
#'
#' @param ... `mm_beta` objects from [fit_beta()] (or a single list of them).
#' @param method `"stouffer"` (default) or `"fisher"`.
#' @return A data.frame with one row per term: the number of studies `k`, the combined
#'   `p_value`, the mean `R2`, and the `method`.
#' @family pooling
#' @export
pool_beta <- function(..., method = c("stouffer", "fisher")) {
  method <- match.arg(method)
  xs <- list(...)
  if (length(xs) == 1 && is.list(xs[[1]]) && !inherits(xs[[1]], "data.frame")) xs <- xs[[1]]
  .check_list_of(xs, "mm_beta", "fit_beta()")
  d <- do.call(rbind, lapply(xs, as.data.frame))
  if (length(unique(d$distance)) > 1) stop("All studies must use the same distance.", call. = FALSE)
  if (anyDuplicated(paste(d$study, d$term))) {
    i <- anyDuplicated(paste(d$study, d$term)); stop("Study '", d$study[i], "' has term '", d$term[i], "' more than once.", call. = FALSE)
  }
  rows <- lapply(unique(d$term), function(tm) {
    x <- d[d$term == tm, ]; k <- nrow(x)
    p <- pmin(pmax(x$p_value, 1e-12), 1 - 1e-12)
    pc <- if (method == "fisher") stats::pchisq(-2 * sum(log(p)), 2 * k, lower.tail = FALSE) else {
      w <- sqrt(x$n); stats::pnorm(sum(w * stats::qnorm(p, lower.tail = FALSE)) / sqrt(sum(w^2)), lower.tail = FALSE)
    }
    data.frame(term = tm, k = k, p_value = pc, R2 = stats::weighted.mean(x$R2, x$n), method = method, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}
