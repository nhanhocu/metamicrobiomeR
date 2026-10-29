# Meta-analysis of per-study effects (see dev/SPEC.md)

#' Pool effects across studies
#'
#' Random-effects meta-analysis of the per-study estimates in an effect table,
#' separately for every feature, model part and term, using `metafor::rma.uni()`.
#'
#' A feature is pooled only if it was estimated in at least `min_studies` studies
#' and in at least `min_prop` of all studies that report that part and term.
#' Studies with a missing estimate or a non-positive standard error are dropped.
#' The Hartung-Knapp tests (`"adhoc"`, `"knha"`) use a t distribution with `k - 1` degrees of
#' freedom, so with only two studies they have one degree of freedom and even a clear effect
#' will not be significant; that is the honest answer with so little replication, and
#' `test = "z"` (which ignores the uncertainty of the between-study variance) is the optimistic one.
#' Pooling is only meaningful when the studies report the same effect scale;
#' [bind_effects()] enforces this.
#'
#' @param effects An `mm_effects` object with several studies (see [bind_effects()]).
#' @param term,part Optional character vectors restricting the terms and model parts.
#' @param method Estimator of the between-study variance, passed to
#'   `metafor::rma.uni()`: "REML" (default), "DL" (DerSimonian-Laird, as in 1.x),
#'   "PM", "EB", "ML", "HE" or "SJ".
#' @param test Test and confidence interval for the pooled effect, passed to
#'   `metafor::rma.uni()`: `"adhoc"` (default) is the Hartung-Knapp adjustment
#'   with the ad hoc correction that keeps the interval from becoming narrower than
#'   the normal-theory one; `"knha"` is plain Hartung-Knapp, which can give
#'   spuriously small p-values when only a few studies agree closely (see
#'   `dev/VALIDATION.md`); `"z"` is the normal-theory test of 1.x. Use
#'   `method = "DL", test = "z"` to reproduce 1.x.
#' @param min_studies,min_prop Minimum number of studies, and minimum proportion
#'   of the studies reporting the part and term (1.x used 0.5).
#' @param p_adjust Method for [stats::p.adjust()], or `"none"`.
#' @param adjust_by Columns defining the family of tests adjusted together: any of
#'   `"part"`, `"term"`, `"rank"`. The default adjusts within each rank, part and
#'   term; `c("part", "term")` mixes all ranks, as 1.x did.
#' @return An object of class `mm_pooled` (a data.frame), one row per feature,
#'   part and term, with the number of studies `k`, the pooled `estimate`, `se`,
#'   confidence limits, test `statistic` and `p_value`, the heterogeneity measures
#'   `tau2`, `I2`, `Q` and `Q_p`, the prediction interval (`pi_lo`, `pi_hi`; `NA`
#'   unless there are at least 3 studies), the share of the total random-effects
#'   weight carried by the most influential study (`max_weight`, between 0 and 1) and
#'   that study (`top_study`), and the adjusted p-value `p_adjust`. A `max_weight`
#'   close to 1 means one study determines the pooled result; use [leave_one_out()]
#'   to see it without that study.
#' @export
pool_effects <- function(effects, term = NULL, part = NULL,
                         method = c("REML", "DL", "PM", "EB", "ML", "HE", "SJ"),
                         test = c("adhoc", "knha", "z"),
                         min_studies = 2, min_prop = 0.5, p_adjust = "fdr",
                         adjust_by = c("part", "term", "rank")) {
  method <- match.arg(method)
  test <- match.arg(test)
  .check_class(effects, "mm_effects", "effects", "fit_study() or bind_effects()")
  check_effects(effects)
  if (!is.null(term)) effects <- effects[effects$term %in% term, , drop = FALSE]
  if (!is.null(part)) effects <- effects[effects$part %in% part, , drop = FALSE]
  if (!nrow(effects)) stop("No effects left after selecting term and part.", call. = FALSE)
  adjust_by <- match.arg(adjust_by, c("part", "term", "rank"), several.ok = TRUE)
  effects <- as.data.frame(effects)
  bad <- is.na(effects$estimate) | is.na(effects$se) | effects$se <= 0
  if (any(bad)) {
    warning(sum(bad), " estimate(s) with a missing value or non-positive standard error were dropped.", call. = FALSE)
    effects <- effects[!bad, , drop = FALSE]
  }
  if (!nrow(effects)) stop("No usable estimates to pool.", call. = FALSE)

  effects$.group <- paste(effects$part, effects$term, sep = "\r")
  n_studies <- tapply(effects$study, effects$.group, function(s) length(unique(s)))
  effects$.key <- paste(effects$.group, effects$feature, sep = "\r")
  keys <- unique(effects$.key)

  rows <- lapply(keys, function(key) {
    d <- effects[effects$.key == key, , drop = FALSE]
    k <- nrow(d)
    if (k < min_studies || k < min_prop * n_studies[[d$.group[1]]]) return(NULL)
    fit <- tryCatch(suppressWarnings(metafor::rma.uni(yi = d$estimate, sei = d$se, method = method, test = test)),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    pi <- tryCatch(suppressWarnings(stats::predict(fit)), error = function(e) NULL)
    rk <- d$rank[!is.na(d$rank)]
    ws <- 1 / (d$se^2 + fit$tau2); share <- ws / sum(ws)
    data.frame(feature = d$feature[1], rank = if (length(rk)) rk[1] else NA_character_,
               part = d$part[1], term = d$term[1], k = k,
               estimate = as.numeric(fit$beta), se = fit$se, ci_lo = fit$ci.lb, ci_hi = fit$ci.ub,
               statistic = fit$zval, p_value = fit$pval, tau2 = fit$tau2, I2 = fit$I2,
               Q = fit$QE, Q_p = fit$QEp,
               pi_lo = if (k >= 3 && !is.null(pi)) pi$pi.lb else NA_real_,
               pi_hi = if (k >= 3 && !is.null(pi)) pi$pi.ub else NA_real_,
               max_weight = max(share), top_study = as.character(d$study[which.max(share)]),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) stop("No feature met the requirements for pooling (min_studies, min_prop).", call. = FALSE)
  fam <- do.call(paste, c(lapply(adjust_by, function(a) out[[a]]), sep = "\r"))
  out$p_adjust <- NA_real_
  if (p_adjust == "none") out$p_adjust <- out$p_value else {
    for (g in unique(fam)) out$p_adjust[fam == g] <- stats::p.adjust(out$p_value[fam == g], method = p_adjust)
  }
  out <- out[order(out$part, out$term, out$p_value), ]
  rownames(out) <- NULL
  attr(out, "settings") <- list(method = method, test = test, min_studies = min_studies, min_prop = min_prop,
                                p_adjust = p_adjust, adjust_by = adjust_by)
  class(out) <- c("mm_pooled", "data.frame")
  out
}

#' Pooled effects with each study left out in turn
#'
#' Repeats [pool_effects()] once per study, omitting that study each time, to show
#' how much a result depends on any single study.
#'
#' @param effects An `mm_effects` object with at least two studies.
#' @param feature Optional character vector restricting the features.
#' @param term,part Optional restrictions, as in [pool_effects()].
#' @param ... Further arguments passed to [pool_effects()] (`method`, `test`, ...).
#'   Adjusted p-values are computed among the features pooled in each run.
#' @return A data.frame with the columns of [pool_effects()] plus `omitted`, the
#'   study left out; features that can no longer be pooled without that study are
#'   absent from its rows.
#' @export
leave_one_out <- function(effects, feature = NULL, term = NULL, part = NULL, ...) {
  .check_class(effects, "mm_effects", "effects", "fit_study() or bind_effects()")
  if (!is.null(feature)) effects <- effects[effects$feature %in% feature, , drop = FALSE]
  studies <- unique(as.character(effects$study))
  if (length(studies) < 2) stop("Leave-one-out needs at least two studies.", call. = FALSE)
  res <- lapply(studies, function(s) {
    p <- tryCatch(suppressWarnings(pool_effects(effects[effects$study != s, , drop = FALSE],
                                                term = term, part = part, ...)), error = function(e) NULL)
    if (is.null(p)) return(NULL)
    p <- as.data.frame(p); p$omitted <- s
    p
  })
  out <- do.call(rbind, res)
  if (is.null(out)) stop("No feature could be pooled after omitting any single study.", call. = FALSE)
  rownames(out) <- NULL
  out
}
