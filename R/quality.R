# Quality diagnostics for per-study effects, and sensitivity of the pooled result

#' Diagnose the per-study estimates before pooling
#'
#' Pooling weights studies by their precision, so an implausibly small standard error
#' gives one study almost all the weight. This compares each study with the others, for each
#' part and term: a study's estimate is flagged when its standard error is below `se_ratio`
#' times the typical small standard error of the other studies (their median 5th percentile).
#' Estimates from fits that did not converge, and (optionally) very large estimates, are flagged
#' too. The comparison needs at least two studies.
#'
#' @param effects An `mm_effects` object with several studies.
#' @param se_ratio Flag standard errors below this fraction of the other studies' typical
#'   small standard error.
#' @param max_abs_estimate Optionally flag estimates with a larger absolute value (only
#'   meaningful on scales with a natural range, such as logits).
#' @param quantile The quantile used as a study's "small standard error" (default 5 per cent).
#' @return An object of class `mm_quality`: `studies` (per study, part and term: the number of
#'   estimates, the proportion converged, the minimum, the chosen quantile and the median of
#'   the standard errors, and how many were flagged) and `rows` (the effect rows with the flag
#'   columns `flag_se`, `flag_converged`, `flag_estimate` and `flagged`).
#' @export
check_effects_quality <- function(effects, se_ratio = 0.1, max_abs_estimate = NULL, quantile = 0.05) {
  .check_class(effects, "mm_effects", "effects", "fit_study() or bind_effects()")
  check_effects(effects)
  d <- as.data.frame(effects)
  d$flag_se <- FALSE; d$flag_converged <- !d$converged
  d$flag_estimate <- if (is.null(max_abs_estimate)) FALSE else !is.na(d$estimate) & abs(d$estimate) > max_abs_estimate
  grp <- paste(d$part, d$term, sep = "\r")
  summ <- list()
  for (g in unique(grp)) {
    ii <- which(grp == g)
    q <- tapply(d$se[ii], d$study[ii], function(s) stats::quantile(s[!is.na(s) & s > 0], quantile, names = FALSE))
    for (s in names(q)) {
      other <- setdiff(names(q), s)
      ref <- if (length(other)) stats::median(q[other], na.rm = TRUE) else NA_real_
      jj <- ii[d$study[ii] == s]
      if (!is.na(ref)) d$flag_se[jj] <- !is.na(d$se[jj]) & d$se[jj] > 0 & d$se[jj] < se_ratio * ref
      se <- d$se[jj]
      summ[[length(summ) + 1]] <- data.frame(study = s, part = d$part[jj][1], term = d$term[jj][1], n = length(jj),
        converged = mean(d$converged[jj]), se_min = suppressWarnings(min(se, na.rm = TRUE)), se_q = unname(q[[s]]),
        se_median = stats::median(se, na.rm = TRUE), reference_se_q = ref, flagged_se = sum(d$flag_se[jj]),
        flagged_estimate = sum(d$flag_estimate[jj]), stringsAsFactors = FALSE)
    }
  }
  d$flagged <- d$flag_se | d$flag_converged | d$flag_estimate
  structure(list(studies = do.call(rbind, summ), rows = d, settings = list(se_ratio = se_ratio, max_abs_estimate = max_abs_estimate, quantile = quantile)),
            class = "mm_quality")
}

#' @export
print.mm_quality <- function(x, ...) {
  s <- x$studies
  cat("<mm_quality> ", nrow(x$rows), " estimates; flagged: ", sum(x$rows$flagged), " (standard error ", sum(x$rows$flag_se),
      ", not converged ", sum(x$rows$flag_converged), ", large estimate ", sum(x$rows$flag_estimate), ")\n", sep = "")
  bad <- s[s$flagged_se > 0 | s$flagged_estimate > 0, , drop = FALSE]
  for (i in seq_len(nrow(bad))) cat(sprintf("  study %s, %s/%s: %d implausibly small standard errors (min %.3g; other studies' typical small value %.3g)\n",
                                            bad$study[i], bad$part[i], bad$term[i], bad$flagged_se[i], bad$se_min[i], bad$reference_se_q[i]))
  invisible(x)
}

#' Drop doubtful estimates before pooling
#'
#' Removes estimates flagged by [check_effects_quality()]: from fits that did not converge
#' and/or with an implausibly small standard error (or very large estimate).
#'
#' @inheritParams check_effects_quality
#' @param drop_nonconverged,drop_implausible_se Which flags to act on.
#' @return The filtered `mm_effects` object, with the removed rows and the reason in the
#'   attribute `dropped`.
#' @export
filter_effects <- function(effects, drop_nonconverged = TRUE, drop_implausible_se = TRUE, se_ratio = 0.1, max_abs_estimate = NULL) {
  q <- check_effects_quality(effects, se_ratio = se_ratio, max_abs_estimate = max_abs_estimate)
  r <- q$rows
  drop <- (drop_nonconverged & r$flag_converged) | (drop_implausible_se & r$flag_se) | (!is.null(max_abs_estimate) & r$flag_estimate)
  out <- as.data.frame(effects)[!drop, , drop = FALSE]
  class(out) <- class(effects); rownames(out) <- NULL
  dropped <- r[drop, , drop = FALSE]
  dropped$reason <- apply(dropped[, c("flag_converged", "flag_se", "flag_estimate")], 1, function(f) paste(c("not converged", "implausible standard error", "large estimate")[f], collapse = ", "))
  attr(out, "dropped") <- dropped[, c("study", "feature", "part", "term", "estimate", "se", "reason")]
  for (a in c("fit_info", "fit_settings")) if (!is.null(attr(effects, a))) attr(out, a) <- attr(effects, a)
  if (nrow(dropped)) message(nrow(dropped), " of ", nrow(r), " estimates dropped before pooling (see attr(, \"dropped\")).")
  out
}

#' Compare pooling settings
#'
#' Pools the same effects under several settings and shows which features are significant
#' under all of them. Results that depend on the estimator or on the test (for example a
#' feature significant with the normal test but not with the Hartung-Knapp test) deserve
#' caution, especially with few studies.
#'
#' @param effects An `mm_effects` object with several studies.
#' @param settings A named list of argument lists for [pool_effects()]. The default compares
#'   the 1.x calculation (`DL_z`), REML with the normal test, REML with plain Hartung-Knapp and
#'   REML with the ad hoc Hartung-Knapp test.
#' @param term,part Optional restrictions, as in [pool_effects()].
#' @param alpha Significance level applied to the adjusted p-value (or to the p-value if
#'   `significance = "p_value"`).
#' @param significance Which p-value defines "significant".
#' @param ... Further arguments passed to [pool_effects()] in every setting (for example
#'   `min_prop` or `adjust_by`); the entries of `settings` take precedence.
#' @return An object of class `mm_poolcompare`: `long` (all pooled rows with the `setting`),
#'   `wide` (one row per feature, part and term with the p-values and the significance under
#'   each setting), `n_significant` (per setting) and `robust` (features significant under all).
#' @export
pool_sensitivity <- function(effects, settings = NULL, term = NULL, part = NULL, alpha = 0.05, significance = c("p_adjust", "p_value"), ...) {
  significance <- match.arg(significance)
  if (is.null(settings)) settings <- list(DL_z = list(method = "DL", test = "z"), REML_z = list(method = "REML", test = "z"),
                                           REML_knha = list(method = "REML", test = "knha"), REML_adhoc = list(method = "REML", test = "adhoc"))
  if (is.null(names(settings)) || any(!nzchar(names(settings)))) stop("'settings' must be a named list.", call. = FALSE)
  long <- do.call(rbind, lapply(names(settings), function(nm) {
    args <- utils::modifyList(c(list(effects = effects, term = term, part = part), list(...)), settings[[nm]])
    p <- do.call(pool_effects, args)
    p <- as.data.frame(p); p$setting <- nm; p
  }))
  long$significant <- !is.na(long[[significance]]) & long[[significance]] <= alpha
  key <- paste(long$part, long$term, long$feature, sep = "\r")
  keys <- unique(key)
  wide <- data.frame(feature = long$feature[match(keys, key)], part = long$part[match(keys, key)], term = long$term[match(keys, key)],
                     k = long$k[match(keys, key)], stringsAsFactors = FALSE)
  for (nm in names(settings)) {
    s <- long[long$setting == nm, ]; m <- match(keys, paste(s$part, s$term, s$feature, sep = "\r"))
    wide[[paste0("estimate_", nm)]] <- s$estimate[m]; wide[[paste0("p_", nm)]] <- s[[significance]][m]; wide[[paste0("sig_", nm)]] <- s$significant[m]
  }
  sigcols <- paste0("sig_", names(settings))
  wide$n_settings_significant <- rowSums(wide[, sigcols, drop = FALSE], na.rm = TRUE)
  wide$robust <- wide$n_settings_significant == length(settings)
  wide <- wide[order(-wide$n_settings_significant, wide$feature), ]; rownames(wide) <- NULL
  structure(list(long = long, wide = wide, n_significant = stats::setNames(colSums(wide[, sigcols, drop = FALSE], na.rm = TRUE), names(settings)),
                 robust = wide$feature[wide$robust], alpha = alpha, significance = significance), class = "mm_poolcompare")
}

#' @export
print.mm_poolcompare <- function(x, ...) {
  cat("<mm_poolcompare> significant (", x$significance, " <= ", x$alpha, ") per setting:\n", sep = "")
  print(x$n_significant)
  cat("significant under all", length(x$n_significant), "settings:", length(x$robust), "of", nrow(x$wide), "features\n")
  invisible(x)
}

#' Combine per-study p-values with a direction
#'
#' An alternative to pooling effect sizes that needs only each study's p-value and the sign of
#' its estimate, so it also works when studies were analysed with different engines or tools
#' whose estimates are not on a common scale. Two-sided p-values are turned into signed z
#' scores and combined with weighted Stouffer's method (weights proportional to the square root
#' of the sample size) or Fisher's method (applied to the one-sided p-values in each direction).
#' It reports the number of studies with positive and negative estimates and a sign test.
#' It gives evidence of an effect, not an estimate of its size, and does not model heterogeneity.
#' Like pooling, it is not protected against degenerate fits: a single study with an extreme
#' p-value (for example from a fit with an implausibly small standard error) can dominate the
#' combination, so run [check_effects_quality()] and [filter_effects()] first.
#'
#' @param effects An `mm_effects` object (studies need not share an engine or scale). The
#'   `p_value` column is used; where it is missing, a normal-approximation p-value is computed
#'   from the estimate and standard error.
#' @param method `"stouffer"` (default) or `"fisher"`.
#' @param weights `"sqrt_n"` (default, Stouffer only) or `"none"`.
#' @param term,part Optional restrictions.
#' @param min_studies Minimum number of studies per feature.
#' @param p_adjust Method for [stats::p.adjust()], or `"none"`.
#' @param adjust_by Columns defining the family adjusted together (`"part"`, `"term"`, `"rank"`).
#' @return A data.frame of class `mm_combined` with `k`, `n_pos`, `n_neg`, `direction`
#'   (`"up"`, `"down"` or `"mixed"`), `sign_test_p`, the combined `statistic` (Stouffer z) and
#'   `p_value`, and `p_adjust`.
#' @export
combine_pvalues <- function(effects, method = c("stouffer", "fisher"), weights = c("sqrt_n", "none"), term = NULL, part = NULL,
                            min_studies = 2, p_adjust = "fdr", adjust_by = c("part", "term", "rank")) {
  method <- match.arg(method); weights <- match.arg(weights); adjust_by <- match.arg(adjust_by, c("part", "term", "rank"), several.ok = TRUE)
  .check_class(effects, "mm_effects", "effects", "fit_study() or bind_effects()")
  check_effects(effects)
  if (!is.null(term)) effects <- effects[effects$term %in% term, , drop = FALSE]
  if (!is.null(part)) effects <- effects[effects$part %in% part, , drop = FALSE]
  d <- as.data.frame(effects)
  miss <- is.na(d$p_value) & !is.na(d$estimate) & !is.na(d$se) & d$se > 0
  d$p_value[miss] <- 2 * stats::pnorm(-abs(d$estimate[miss] / d$se[miss]))
  d <- d[!is.na(d$p_value) & !is.na(d$estimate), , drop = FALSE]
  if (!nrow(d)) stop("No estimates with a p-value.", call. = FALSE)
  key <- paste(d$part, d$term, d$feature, sep = "\r")
  rows <- lapply(unique(key), function(k) {
    x <- d[key == k, , drop = FALSE]; K <- nrow(x)
    if (K < min_studies) return(NULL)
    p <- pmin(pmax(x$p_value, 1e-300), 1); sg <- sign(x$estimate)
    z <- NA_real_
    if (method == "stouffer") {
      w <- if (weights == "sqrt_n") sqrt(ifelse(is.na(x$n) | x$n < 1, 1, x$n)) else rep(1, K)
      z <- sum(w * sg * stats::qnorm(p / 2, lower.tail = FALSE)) / sqrt(sum(w^2))
      pc <- 2 * stats::pnorm(-abs(z))
    } else {
      pp <- ifelse(sg >= 0, p / 2, 1 - p / 2)
      pc <- min(1, 2 * min(stats::pchisq(-2 * sum(log(pp)), 2 * K, lower.tail = FALSE), stats::pchisq(-2 * sum(log(1 - pp)), 2 * K, lower.tail = FALSE)))
    }
    np <- sum(sg > 0); nn <- sum(sg < 0)
    rk <- x$rank[!is.na(x$rank)]
    data.frame(feature = x$feature[1], rank = if (length(rk)) rk[1] else NA_character_, part = x$part[1], term = x$term[1], k = K,
               n_pos = np, n_neg = nn, direction = if (nn == 0 && np > 0) "up" else if (np == 0 && nn > 0) "down" else "mixed",
               sign_test_p = if (np + nn > 0) stats::binom.test(np, np + nn, 0.5)$p.value else NA_real_,
               statistic = z, p_value = pc, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) stop("No feature has estimates from at least ", min_studies, " studies.", call. = FALSE)
  fam <- do.call(paste, c(lapply(adjust_by, function(a) out[[a]]), sep = "\r")); out$p_adjust <- NA_real_
  for (g in unique(fam)) out$p_adjust[fam == g] <- if (p_adjust == "none") out$p_value[fam == g] else stats::p.adjust(out$p_value[fam == g], p_adjust)
  out <- out[order(out$part, out$term, out$p_value), ]; rownames(out) <- NULL
  attr(out, "settings") <- list(method = method, weights = weights, p_adjust = p_adjust, adjust_by = adjust_by)
  class(out) <- c("mm_combined", "data.frame")
  out
}
