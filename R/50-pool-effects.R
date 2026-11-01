# Pooling effect tables across studies, and leaving studies out

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
#'   `metafor::rma.uni()`: "REML" (default), "DL" (DerSimonian-Laird),
#'   "PM", "EB", "ML", "HE" or "SJ".
#' @param test Test and confidence interval for the pooled effect, passed to
#'   `metafor::rma.uni()`: `"adhoc"` (default) is the Hartung-Knapp adjustment
#'   with the ad hoc correction that keeps the interval from becoming narrower than
#'   the normal-theory one; `"knha"` is plain Hartung-Knapp, which can give
#'   spuriously small p-values when only a few studies agree closely; `"z"` is the
#'   normal-theory test, which ignores the uncertainty of the between-study variance and is
#'   anti-conservative with few studies.
#' @param min_studies,min_prop Minimum number of studies, and minimum proportion
#'   of the studies reporting the part and term.
#' @param p_adjust Method for [stats::p.adjust()], or `"none"`.
#' @param adjust_by Columns defining the family of tests adjusted together: any of
#'   `"part"`, `"term"`, `"rank"`. The default adjusts within each rank, part and
#'   term; `c("part", "term")` mixes all ranks.
#' @return An object of class `mm_pooled` (a data.frame), one row per feature,
#'   part and term, with the number of studies `k`, the pooled `estimate`, `se`,
#'   confidence limits, test `statistic` and `p_value`, the heterogeneity measures
#'   `tau2`, `I2`, `Q` and `Q_p`, the prediction interval (`pi_lo`, `pi_hi`; `NA`
#'   unless there are at least 3 studies), the share of the total random-effects
#'   weight carried by the most influential study (`max_weight`, between 0 and 1) and
#'   that study (`top_study`), and the adjusted p-value `p_adjust`. A `max_weight`
#'   close to 1 means one study determines the pooled result; use [leave_one_out()]
#'   to see it without that study.
#' @family pooling
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
  if (!nrow(effects)) stop("No effects left after selecting term", if (!is.null(term)) paste0(" '", paste(term, collapse = "', '"), "'"), " and part",
                           if (!is.null(part)) paste0(" '", paste(part, collapse = "', '"), "'"), ".", call. = FALSE)
  .check_comparable(effects)
  adjust_by <- match.arg(adjust_by, c("part", "term", "rank"), several.ok = TRUE)
  effects <- as.data.frame(effects)
  miss <- is.na(effects$estimate) | is.na(effects$se); bad <- miss | effects$se <= 0
  if (any(bad)) {
    why <- c(if (any(miss)) paste(sum(miss), "with a missing estimate or standard error"),
             if (any(bad & !miss)) paste(sum(bad & !miss), "with a standard error that is not positive"))
    warning(.plural(sum(bad), "estimate"), " not used: ", paste(why, collapse = ", "), ".", call. = FALSE)
    effects <- effects[!bad, , drop = FALSE]
  }
  if (!nrow(effects)) stop("No usable estimates to pool: every estimate has a missing value or a standard error that is missing or not positive.", call. = FALSE)

  effects$.group <- paste(effects$part, effects$term, sep = "\r")
  n_studies <- tapply(effects$study, effects$.group, function(s) length(unique(s)))
  effects$.key <- paste(effects$.group, effects$feature, sep = "\r")
  keys <- unique(effects$.key)

  est <- effects$estimate; sev <- effects$se; stu <- as.character(effects$study); rk_all <- effects$rank
  grp_idx <- split(seq_len(nrow(effects)), factor(effects$.key, levels = keys))      # one pass, not one scan per feature
  fast <- method %in% c("REML", "DL", "PM")
  pool_one <- function(ii) {
    k <- length(ii)
    if (k < min_studies || k < min_prop * n_studies[[effects$.group[ii[1]]]]) return(NULL)
    y <- est[ii]; v <- sev[ii]^2
    r <- if (fast) .rma_fast(y, v, method = method, test = test) else .rma_metafor(y, v, method, test)
    if (is.null(r)) return(NULL)
    ws <- 1 / (v + r$tau2); share <- ws / sum(ws)
    rk <- rk_all[ii]; rk <- rk[!is.na(rk)]
    list(feature = effects$feature[ii[1]], rank = if (length(rk)) rk[1] else NA_character_, part = effects$part[ii[1]], term = effects$term[ii[1]],
         k = k, estimate = r$beta[1], se = r$se[1], ci_lo = r$ci_lb[1], ci_hi = r$ci_ub[1], statistic = r$stat[1], p_value = r$pval[1],
         tau2 = r$tau2, I2 = r$I2, Q = r$QE, Q_p = r$QEp,
         pi_lo = if (k >= 3 && !is.null(r$pi_lb)) r$pi_lb else NA_real_, pi_hi = if (k >= 3 && !is.null(r$pi_ub)) r$pi_ub else NA_real_,
         max_weight = max(share), top_study = stu[ii][which.max(share)])
  }
  rows <- lapply(grp_idx, pool_one)
  rows <- rows[!vapply(rows, is.null, logical(1))]
  out <- if (length(rows)) {
    cols <- names(rows[[1]])
    as.data.frame(lapply(stats::setNames(cols, cols), function(cn) unlist(lapply(rows, `[[`, cn), use.names = FALSE)), stringsAsFactors = FALSE)
  } else NULL
  if (is.null(out)) {
    kmax <- max(table(effects$.key))
    stop("No feature was estimated in enough studies to pool: it needs at least ", min_studies, " studies and ", .pct(min_prop),
         " of the studies reporting the term, but the most any feature has is ", kmax, ". Lower 'min_studies' or 'min_prop'.", call. = FALSE)
  }
  n_total <- length(keys)
  if (nrow(out) < n_total) {
    .inform(.plural(n_total - nrow(out), "feature"), " of ", n_total, " not pooled (estimated in fewer than ", min_studies, " studies",
            if (min_prop > 0) paste0(" or in fewer than ", .pct(min_prop), " of the studies reporting the term"), ").")
  }
  if (test != "z" && any(out$k == 2)) {
    .inform(.plural(sum(out$k == 2), "feature"), if (sum(out$k == 2) == 1) " was" else " were", " pooled from only 2 studies: the Hartung-Knapp test then has 1 degree of freedom and is very conservative (see ?pool_effects).",
            once = "k2-hk")
  }
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
#' @family pooling
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
  if (is.null(out)) stop("No feature could be pooled after omitting any single study; each needs at least 'min_studies' studies left (default 2).", call. = FALSE)
  rownames(out) <- NULL
  out
}

#' @export
print.mm_pooled <- function(x, n = 6, ...) {
  d <- as.data.frame(x); s <- attr(x, "settings")
  if (!all(c("feature", "part", "term", "k", "estimate", "p_value", "I2") %in% colnames(d)) || !nrow(d)) { print(d, ...); return(invisible(x)) }
  cat("<mm_pooled> ", .plural(nrow(d), "feature"), " pooled",
      if (!is.null(s)) paste0(" (", s$method, ", ", s$test, " test)"), "; ", if (min(d$k) == max(d$k)) max(d$k) else paste(min(d$k), "to", max(d$k)), " studies each\n", sep = "")
  sig <- if ("p_adjust" %in% colnames(d)) sum(d$p_adjust <= 0.05, na.rm = TRUE) else NA
  one_pt <- length(unique(d$part)) == 1 && length(unique(d$term)) == 1
  cat("  ", if (one_pt) paste0("term ", d$term[1], " (", d$part[1], ") | "), "adjusted p <= 0.05: ", sig, " | median I2: ", format(round(stats::median(d$I2, na.rm = TRUE)), trim = TRUE), "%",
      if ("max_weight" %in% colnames(d) && any(d$max_weight > 0.9, na.rm = TRUE)) paste0(" | ", sum(d$max_weight > 0.9, na.rm = TRUE), " dominated by one study (max_weight > 0.9)"), "\n", sep = "")
  cols <- c("feature", if (!one_pt) c("part", "term"), "k", "estimate", "ci_lo", "ci_hi", "p_value", "p_adjust", "I2")
  show <- d[seq_len(min(n, nrow(d))), intersect(cols, colnames(d)), drop = FALSE]
  print(show, row.names = FALSE, digits = 3)
  if (nrow(d) > n) cat("# ... ", nrow(d) - n, " more rows\n", sep = "")
  invisible(x)
}

# metafor, for the estimators the fast core does not cover (EB, ML, HE, SJ)
.rma_metafor <- function(y, v, method, test) {
  .need("metafor", "the estimators EB, ML, HE and SJ (REML, DL and PM need no extra package)")
  fit <- tryCatch(suppressWarnings(metafor::rma.uni(yi = y, vi = v, method = method, test = test)), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  pr <- if (length(y) >= 3) tryCatch(suppressWarnings(stats::predict(fit)), error = function(e) NULL) else NULL
  list(beta = as.numeric(fit$beta), se = fit$se, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub, stat = fit$zval, pval = fit$pval, tau2 = fit$tau2,
       I2 = fit$I2, QE = fit$QE, QEp = fit$QEp, pi_lb = if (!is.null(pr)) pr$pi.lb, pi_ub = if (!is.null(pr)) pr$pi.ub)
}
