# Beta diversity: PERMANOVA per study and combination across studies

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
#' @export
fit_beta <- function(study, formula, distance = "bray", permutations = 999, rank = NULL, seed = NULL,
                     longitudinal = NULL) {
  .need("vegan", "fit_beta()"); .need("permute", "fit_beta()")
  if (!inherits(study, "mm_study")) stop("'study' must be an mm_study object.", call. = FALSE)
  if (!inherits(formula, "formula") || length(formula) != 2) stop("'formula' must be one-sided, e.g. ~ group.", call. = FALSE)
  info <- mm_study_info(study); sdat <- mm_samples(study)
  vars <- all.vars(formula); miss <- setdiff(vars, colnames(sdat))
  if (length(miss)) stop("Variable(s) not in the sample data: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  ra <- mm_relabund(study); feats <- mm_features(study)
  if (!is.null(rank)) {
    if (!"rank" %in% colnames(feats)) stop("The feature data has no 'rank' column.", call. = FALSE)
    ra <- ra[feats$rank %in% rank, , drop = FALSE]
    if (!nrow(ra)) stop("No features with rank: ", paste(rank, collapse = ", "), ".", call. = FALSE)
  }
  d0 <- sdat[, vars, drop = FALSE]
  ok <- stats::complete.cases(d0) & colSums(is.na(ra)) == 0 & colSums(ra) > 0
  if (sum(ok) < 4) stop("Too few complete samples.", call. = FALSE)
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
    if (length(varying)) stop("Term(s) vary within subjects: ", paste(varying, collapse = ", "),
                              ". Permuting whole subjects can only test between-subject terms. ",
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
#' @export
pool_beta <- function(..., method = c("stouffer", "fisher")) {
  method <- match.arg(method)
  xs <- list(...)
  if (length(xs) == 1 && is.list(xs[[1]]) && !inherits(xs[[1]], "data.frame")) xs <- xs[[1]]
  if (!length(xs) || !all(vapply(xs, inherits, logical(1), "mm_beta"))) stop("Supply mm_beta objects (see fit_beta()).", call. = FALSE)
  d <- do.call(rbind, lapply(xs, as.data.frame))
  if (length(unique(d$distance)) > 1) stop("All studies must use the same distance.", call. = FALSE)
  if (anyDuplicated(paste(d$study, d$term))) stop("Duplicate study and term.", call. = FALSE)
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
