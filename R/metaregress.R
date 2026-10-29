# Meta-regression: do study-level characteristics explain differences between studies?

#' Meta-regression on study-level moderators
#'
#' For every feature, model part and term, regresses the per-study estimates on
#' study-level characteristics (for example sequencing region, platform or
#' database) with `metafor::rma.uni()`. This asks whether the effect differs
#' between kinds of studies. With few studies the power is low, and a feature
#' needs more studies than there are model coefficients.
#'
#' @param effects An `mm_effects` object with several studies (see [bind_effects()]).
#' @param study_data A data.frame with a column `study` and one column per moderator
#'   (see [study_table()]).
#' @param moderators One-sided formula of the moderators, e.g. `~ region`.
#' @param term,part Optional restrictions, as in [pool_effects()].
#' @param method,test As in [pool_effects()].
#' @param p_adjust Method for [stats::p.adjust()] applied to the omnibus moderator test
#'   within each rank, part and term, or `"none"`.
#' @return A data.frame of class `mm_metareg`, one row per feature, part, term and
#'   moderator coefficient (the intercept is excluded), with the number of studies `k`,
#'   the coefficient `estimate`, `se`, confidence limits, `statistic` and `p_value`, the
#'   omnibus test of all moderators (`QM`, `QM_p`, adjusted `QM_p_adjust`), the residual
#'   heterogeneity `tau2` and `I2`.
#' @export
meta_regress <- function(effects, study_data, moderators, term = NULL, part = NULL,
                         method = c("REML", "DL", "PM", "EB", "ML", "HE", "SJ"),
                         test = c("adhoc", "knha", "z"), p_adjust = "fdr") {
  method <- match.arg(method); test <- match.arg(test)
  .check_class(effects, "mm_effects", "effects", "fit_study() or bind_effects()")
  check_effects(effects)
  if (!is.data.frame(study_data) || !"study" %in% colnames(study_data)) {
    stop("'study_data' must be a data.frame with a column 'study'.", call. = FALSE)
  }
  if (!inherits(moderators, "formula") || length(moderators) != 2) stop("'moderators' must be a one-sided formula, e.g. ~ region.", call. = FALSE)
  mv <- all.vars(moderators)
  miss <- setdiff(mv, colnames(study_data))
  if (length(miss)) stop("Moderator(s) not in 'study_data': ", paste(miss, collapse = ", "), ".", call. = FALSE)
  if (anyDuplicated(study_data$study)) stop("'study_data' must have one row per study.", call. = FALSE)
  if (!is.null(term)) effects <- effects[effects$term %in% term, , drop = FALSE]
  if (!is.null(part)) effects <- effects[effects$part %in% part, , drop = FALSE]
  .check_comparable(effects)
  effects <- as.data.frame(effects)
  unknown <- setdiff(unique(effects$study), study_data$study)
  if (length(unknown)) stop("No study_data for: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  effects <- effects[!is.na(effects$estimate) & !is.na(effects$se) & effects$se > 0, , drop = FALSE]
  if (!nrow(effects)) stop("No usable estimates.", call. = FALSE)
  effects <- cbind(effects, study_data[match(effects$study, study_data$study), mv, drop = FALSE])
  key <- paste(effects$part, effects$term, effects$feature, sep = "\r")
  rows <- lapply(unique(key), function(k) {
    d <- effects[key == k, , drop = FALSE]
    for (v in mv) if (is.character(d[[v]])) d[[v]] <- factor(d[[v]])
    for (v in mv) if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
    if (any(vapply(mv, function(v) length(unique(d[[v]])) < 2, logical(1)))) return(NULL)   # a moderator without variation
    mm <- stats::model.matrix(moderators, d)
    if (nrow(d) <= ncol(mm) || qr(mm)$rank < ncol(mm)) return(NULL)                          # too few studies or collinear
    fit <- tryCatch(suppressWarnings(metafor::rma.uni(yi = d$estimate, sei = d$se, mods = moderators, data = d,
                                                      method = method, test = test)), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    cf <- seq_along(fit$beta)[-1]
    rk <- d$rank[!is.na(d$rank)]
    data.frame(feature = d$feature[1], rank = if (length(rk)) rk[1] else NA_character_, part = d$part[1], term = d$term[1],
               k = nrow(d), moderator = rownames(fit$beta)[cf], estimate = as.numeric(fit$beta)[cf], se = fit$se[cf],
               ci_lo = fit$ci.lb[cf], ci_hi = fit$ci.ub[cf], statistic = fit$zval[cf], p_value = fit$pval[cf],
               QM = fit$QM, QM_p = fit$QMp, tau2 = fit$tau2, I2 = fit$I2, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  n_skip <- sum(vapply(rows, is.null, logical(1)))
  if (is.null(out)) stop("No feature has enough studies for this meta-regression: it needs more studies than model coefficients and a moderator that varies between them. Check 'study_data' and the formula.", call. = FALSE)
  if (n_skip) .inform(n_skip, " of ", length(rows), " features skipped (fewer studies than model coefficients, or a moderator without variation).")
  one <- !duplicated(paste(out$part, out$term, out$feature, sep = "\r"))
  fam <- paste(out$rank, out$part, out$term, sep = "\r")
  out$QM_p_adjust <- NA_real_
  for (g in unique(fam)) {
    ii <- which(fam == g & one); adj <- if (p_adjust == "none") out$QM_p[ii] else stats::p.adjust(out$QM_p[ii], method = p_adjust)
    for (j in seq_along(ii)) out$QM_p_adjust[fam == g & paste(out$part, out$term, out$feature, sep = "\r") ==
                                                paste(out$part[ii[j]], out$term[ii[j]], out$feature[ii[j]], sep = "\r")] <- adj[j]
  }
  out <- out[order(out$part, out$term, out$QM_p), ]
  rownames(out) <- NULL
  attr(out, "settings") <- list(method = method, test = test, moderators = deparse(moderators), p_adjust = p_adjust)
  class(out) <- c("mm_metareg", "data.frame")
  out
}
