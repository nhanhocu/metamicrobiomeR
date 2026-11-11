# The front door: from a list of studies to pooled results in one call

# Re-signal messages and warnings from one study with the study's name in front
.with_study <- function(id, expr) {
  withCallingHandlers(expr,
    message = function(m) { .inform("[", id, "] ", sub("\n$", "", conditionMessage(m))); invokeRestart("muffleMessage") },    # re-signal first: invokeRestart() does not return
    warning = function(w) { warning("[", id, "] ", conditionMessage(w), call. = FALSE); invokeRestart("muffleWarning") })
}

#' Meta-analysis of several studies in one call
#'
#' Runs the whole workflow with safe defaults: checks that the studies can be pooled
#' ([check_poolable()]), optionally agglomerates them to a rank ([agglomerate()], matching taxa by
#' name across databases), fits the same model to every study ([fit_study()]), checks the estimates
#' ([check_effects_quality()]) and pools them ([pool_effects()]). Messages and warnings are labelled
#' with the study they come from. The steps remain available separately for full control; the
#' returned object holds every intermediate result. `summary()` describes how much the studies differ
#' (heterogeneity), whether conclusions would hold in a new study, and what can explain differences.
#'
#' @param studies A list of `mm_study` objects.
#' @param formula One-sided model formula, e.g. `~ group + age`.
#' @param term Model term(s) to pool (for example `"groupb"`); by default all terms.
#' @param engine A [fit_study()] engine, or `"auto"` (default). With repeated samples (a study with a
#'   subject column, unless the fitting arguments set `longitudinal = FALSE`) it is the centred log-ratio linear
#'   mixed model (`"lm_clr"`), which kept false positives near 5 per cent in every simulated setting with
#'   repeated samples, while the count and beta mixed models exceeded 7.5 per cent with sparse taxa and the
#'   count model often failed to converge. Otherwise it is the count model if every study has counts,
#'   no study has a median library size above one million (as in shotgun profiles) and `glmmTMB` is
#'   installed, and the fast arcsine-square-root linear model if not.
#' @param rank Optional taxonomic rank to analyse. Studies that are not yet agglomerated are
#'   agglomerated to it, matching taxa by harmonised name.
#' @param method,test Pooling settings, see [pool_effects()].
#' @param quality What to do about doubtful estimates (fits that did not converge, and standard errors far
#'   smaller than in the other studies): `"filter"` (default) drops them before pooling and says how many;
#'   `"warn"` keeps them and warns; `"off"` skips the check. Fits that did not converge often have missing or
#'   near-zero standard errors, which can give one study almost all the weight, so they should not be pooled.
#' @param check Run [check_poolable()] first? Problems are reported as warnings, notes as messages.
#' @param fit_args,pool_args Further arguments for [fit_study()] and [pool_effects()] as named lists.
#' @param cores Cores for fitting, see [fit_study()].
#' @param alpha,beta Also analyse alpha diversity ([fit_alpha()]) or beta diversity ([fit_beta()] and
#'   [pool_beta()]) with the same model? Diversity is computed at the analysed level (the `rank`, if given).
#'   Beta diversity needs the `vegan` package and does not allow terms that vary within a subject.
#' @param alpha_args,beta_args Further arguments for [fit_alpha()] and [fit_beta()] as named lists
#'   (for example `list(distance = "aitchison", permutations = 999)`).
#' @return An object of class `mm_meta` with the pooled table (`$pooled`), the per-study effects
#'   (`$effects`), the study table, the checks and the settings. It prints a summary,
#'   `plot()` draws a forest plot (or a heatmap or volcano plot), and `as.data.frame()` gives the pooled table.
#' @examples
#' res <- meta_analyse(example_studies, ~ group + age, term = "groupcase", rank = "genus",
#'                     engine = "lm_asin", fit_args = list(longitudinal = FALSE))
#' res
#' @family pooling
#' @export
meta_analyse <- function(studies, formula, term = NULL, engine = c("auto", .mm_engines), rank = NULL, method = "REML", test = "adhoc",
                         quality = c("filter", "warn", "off"), check = TRUE, fit_args = list(), pool_args = list(),
                         cores = getOption("metamicrobiomeR.cores", 1L), alpha = FALSE, beta = FALSE, alpha_args = list(), beta_args = list()) {
  engine <- match.arg(engine); quality <- match.arg(quality)
  xs <- .as_study_list(list(studies))
  if (length(xs) < 2) stop("Meta-analysis needs at least two studies; got ", length(xs), ". Use fit_study() for a single study.", call. = FALSE)
  ids <- vapply(xs, function(x) x$study$id, character(1))
  if (anyDuplicated(ids)) stop("Study IDs must be unique; '", ids[anyDuplicated(ids)], "' is used twice.", call. = FALSE)
  if (!inherits(formula, "formula") || length(formula) != 2) stop("'formula' must be one-sided, e.g. ~ group + age.", call. = FALSE)
  if (engine == "auto") {
    # repeated samples: as in fit_study(), unless fit_args says otherwise, a study with a subject column is fitted with a subject random intercept
    long <- if (!is.null(fit_args$longitudinal)) isTRUE(fit_args$longitudinal) else any(vapply(xs, function(x) !is.na(mm_study_info(x)$subject), logical(1)))
    all_counts <- all(vapply(xs, function(x) !is.null(mm_counts(x)), logical(1)))
    deep <- all_counts && any(vapply(xs, .large_depth, logical(1)))
    ok <- all_counts && !deep && requireNamespace("glmmTMB", quietly = TRUE)
    engine <- if (long) "lm_clr" else if (ok) "zibb_glmmtmb" else "lm_asin"
    .inform("Using the '", engine, "' engine (", if (long) "repeated samples: in simulations it kept false positives near 5% with sparse taxa, where the count and beta mixed models did not" else if (ok) "every study has counts; about a second per feature and study, engine = \"lm_asin\" is much faster for a first look" else if (!all_counts) "not every study has counts" else if (deep) "a study has a median library size above one million, as in shotgun profiles, where the count model often fails to converge" else "glmmTMB is not installed", ").")
  }
  pc <- NULL
  if (check) {
    pc <- check_poolable(xs, rank = if (!is.null(rank) && all(vapply(xs, function(x) rank %in% colnames(mm_features(x)), logical(1)))) rank)
    w <- pc$issues[pc$issues$level == "warning", , drop = FALSE]; nt <- pc$issues[pc$issues$level == "note", , drop = FALSE]
    if (nrow(w)) warning("The studies may not be safe to pool: ", paste0(w$issue, " (", w$detail, ")", collapse = "; "), ". See check_poolable().", call. = FALSE)
    if (nrow(nt)) .inform("Note: ", paste(nt$issue, collapse = "; "), ". See check_poolable().")
  }
  # small standard errors are handled by the 'quality' step below, so the warnings of fit_study() and pool_effects() about them are not repeated
  quiet_se <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (grepl("implausibly small standard error|standard error below 2% of the median", conditionMessage(w))) invokeRestart("muffleWarning")
  })
  fits <- vector("list", length(xs)); prep <- xs
  for (i in seq_along(xs)) {
    x <- xs[[i]]
    .inform("Study ", i, " of ", length(xs), ": ", ids[i])
    fits[[i]] <- tryCatch(.with_study(ids[i], {
      if (!is.null(rank)) {
        f <- mm_features(x)
        if (!("rank" %in% colnames(f) && any(f$rank %in% rank))) x <- agglomerate(x, rank, key = "name")
      }
      prep[[i]] <- x
      quiet_se(do.call(fit_study, c(list(study = x, formula = formula, engine = engine, rank = rank, cores = cores), fit_args)))
    }), error = function(e) stop("Study '", ids[i], "': ", conditionMessage(e), call. = FALSE))
  }
  bound <- bind_effects(fits)
  q <- NULL
  if (quality != "off") {
    q <- check_effects_quality(bound)
    n_flag <- sum(q$rows$flagged)
    if (n_flag && quality == "warn") {
      by_study <- table(q$rows$study[q$rows$flagged])
      warning(n_flag, " of ", nrow(bound), " estimates look doubtful (", paste0(names(by_study), ": ", by_study, collapse = ", "),
              "): a standard error far smaller than in the other studies or the rest of the same study, or a fit that did not converge. They are kept; use quality = \"filter\" to drop them, ",
              "or see check_effects_quality().", call. = FALSE)
    } else if (n_flag && quality == "filter") bound <- filter_effects(bound)
  }
  pooled <- quiet_se(do.call(pool_effects, c(list(effects = bound, term = term, method = method, test = test), pool_args)))
  per_study <- function(fun, what, args) lapply(seq_along(prep), function(i) tryCatch(.with_study(ids[i], do.call(fun, c(list(study = prep[[i]], formula = formula), args))),
                                                                                   error = function(e) stop("Study '", ids[i], "' (", what, "): ", conditionMessage(e), call. = FALSE)))
  alpha_effects <- alpha_pooled <- beta_pooled <- beta_fits <- NULL
  if (alpha) {
    a_args <- utils::modifyList(if (!is.null(fit_args$longitudinal)) list(longitudinal = fit_args$longitudinal) else list(), alpha_args)
    alpha_effects <- bind_effects(per_study(fit_alpha, "alpha diversity", a_args))
    alpha_pooled <- quiet_se(do.call(pool_effects, c(list(effects = alpha_effects, term = term, method = method, test = test), pool_args[setdiff(names(pool_args), "min_prop")], list(min_prop = 0))))
  }
  if (beta) {
    .need("vegan", "beta diversity (PERMANOVA)")
    beta_fits <- per_study(fit_beta, "beta diversity", utils::modifyList(list(permutations = 499), beta_args))
    beta_pooled <- do.call(pool_beta, beta_fits)
  }
  structure(list(pooled = pooled, effects = bound, alpha = alpha_pooled, alpha_effects = alpha_effects, beta = beta_pooled, beta_fits = beta_fits, studies = study_table(xs), checks = pc, quality = q,
                 settings = list(formula = deparse(formula), engine = engine, rank = rank, term = term, method = method, test = test, quality = quality, alpha = alpha, beta = beta)),
            class = "mm_meta")
}

#' @export
print.mm_meta <- function(x, n = 6, ...) {
  s <- x$settings; p <- as.data.frame(x$pooled)
  cat("<mm_meta> ", .plural(nrow(x$studies), "study", "studies"), " (", .show(x$studies$study, 5), ") | engine ", s$engine, " | model ", s$formula,
      if (!is.null(s$rank)) paste0(" | rank ", s$rank), "\n", sep = "")
  sig <- if ("p_adjust" %in% colnames(p)) sum(p$p_adjust <= 0.05, na.rm = TRUE) else NA
  cat("  ", .plural(nrow(p), "feature"), " pooled (", s$method, ", ", s$test, " test) | adjusted p <= 0.05: ", sig,
      " | median I2: ", format(round(stats::median(p$I2, na.rm = TRUE)), trim = TRUE), "%\n", sep = "")
  dom <- sum(p$max_weight > 0.9, na.rm = TRUE)
  flag <- if (!is.null(x$quality)) sum(x$quality$rows$flagged) else NA
  cat("  checks: ", if (!is.null(x$checks)) paste0(sum(x$checks$issues$level == "warning"), " warning(s), ", sum(x$checks$issues$level == "note"), " note(s)") else "not run",
      " | doubtful estimates: ", if (is.na(flag)) "not checked" else flag, if (dom) paste0(" | ", dom, " dominated by one study") else "", "\n", sep = "")
  top <- p[order(p$p_value), , drop = FALSE]
  one <- length(unique(top$part)) == 1 && length(unique(top$term)) == 1
  cols <- c("feature", if (!one) c("part", "term"), "k", "estimate", "ci_lo", "ci_hi", "p_value", "p_adjust", "I2")
  print(top[seq_len(min(n, nrow(top))), intersect(cols, colnames(top)), drop = FALSE], row.names = FALSE, digits = 3)
  if (nrow(top) > n) cat("# ... ", nrow(top) - n, " more features\n", sep = "")
  if (!is.na(sig) && sig == 0) {
    cat("No feature reaches adjusted p <= 0.05; the table is ordered by p-value. With ", max(p$k), " studies there is little power to pass the correction ",
        "(see pool_sensitivity(x$effects) for how the result depends on the test).\n", sep = "")
  }
  hi <- stats::median(p$I2, na.rm = TRUE)
  if (!is.na(hi) && hi > 50) cat("Heterogeneity is high (median I2 ", round(hi), "%): the studies disagree. summary(x) shows how much and what can explain it.\n", sep = "")
  if (!is.null(x$alpha)) {
    a <- as.data.frame(x$alpha)
    cat("Alpha diversity (pooled standardised difference): ", paste0(a$feature, " ", sprintf("%+.2f", a$estimate), " (p ", format.pval(a$p_value, digits = 2), ")", collapse = "; "), "\n", sep = "")
  }
  if (!is.null(x$beta)) {
    b <- x$beta
    cat("Beta diversity (PERMANOVA, combined across studies): ", paste0(b$term, " R2 ", sprintf("%.2f", b$R2), ", p ", format.pval(b$p_value, digits = 2), collapse = "; "), "\n", sep = "")
  }
  cat("plot(x) draws a forest plot; summary(x) describes the between-study variation; as.data.frame(x) is the pooled table.\n")
  invisible(x)
}

#' @param object An `mm_meta` object.
#' @param top Number of most heterogeneous features to list.
#' @rdname meta_analyse
#' @export
summary.mm_meta <- function(object, top = 5, ...) {
  p <- as.data.frame(object$pooled)
  excl0 <- !is.na(p$ci_lo) & (p$ci_lo > 0 | p$ci_hi < 0)
  pi_ok <- !is.na(p$pi_lo)
  pi0 <- pi_ok & p$pi_lo <= 0 & p$pi_hi >= 0
  ord <- order(-p$I2)
  structure(list(n_features = nrow(p), k = range(p$k), I2 = stats::quantile(p$I2, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE),
                 tau2_median = stats::median(p$tau2, na.rm = TRUE), share_Q = mean(p$Q_p < 0.05, na.rm = TRUE),
                 n_excl0 = sum(excl0), n_not_replicable = sum(excl0 & pi0), n_pi = sum(excl0 & pi_ok),
                 top = p[ord[seq_len(min(top, nrow(p)))], intersect(c("feature", "term", "k", "estimate", "I2", "tau2"), colnames(p)), drop = FALSE],
                 influential = sum(p$max_weight > 0.9, na.rm = TRUE), settings = object$settings),
            class = "summary.mm_meta")
}

#' @export
print.summary.mm_meta <- function(x, ...) {
  cat("<summary of mm_meta> ", .plural(x$n_features, "feature"), ", ", if (x$k[1] == x$k[2]) x$k[1] else paste(x$k, collapse = " to "), " studies each\n\n", sep = "")
  cat("How much do the studies differ?\n")
  cat("  I2 (share of variation beyond sampling error): median ", round(x$I2[2]), "% (middle half ", round(x$I2[1]), "% to ", round(x$I2[3]), "%)\n", sep = "")
  cat("  ", .pct(x$share_Q), " of features have a significant heterogeneity test (Q p < 0.05); median between-study variance tau2 = ", format(signif(x$tau2_median, 2)), "\n", sep = "")
  cat("  Rough convention: below 25% low, 25 to 50% moderate, above 50% substantial. With few studies I2 is imprecise.\n\n")
  cat("Do the conclusions hold in a new study?\n")
  if (x$n_pi > 0) {
    cat("  ", x$n_excl0, " features have a pooled interval that excludes zero; for ", x$n_not_replicable,
        " the prediction interval still includes zero: present on average, but a new study could show no effect.\n", sep = "")
  } else cat("  No feature has a pooled interval that excludes zero, or there are too few studies for prediction intervals.\n")
  if (x$influential > 0) cat("  ", .plural(x$influential, "feature"), " dominated by one study (weight above 90%): see check_effects_quality().\n", sep = "")
  cat("\nMost heterogeneous:\n"); print(x$top, row.names = FALSE, digits = 3)
  cat("\nWhat can explain differences between studies:\n")
  cat("  * Population, protocol or platform: test it with meta_regress() and study characteristics.\n")
  cat("  * One influential study: leave_one_out().\n")
  cat("  * Batch effects that differ between groups within a study (study-wide shifts cancel and do not raise I2):\n")
  cat("    add the batch as a covariate (~ group + batch) in fit_study(), or investigate with other packages and re-run.\n")
  cat("  * Different methods or databases: check_poolable().\n")
  invisible(x)
}

#' @export
as.data.frame.mm_meta <- function(x, ...) as.data.frame(x$pooled, ...)

#' @param x An `mm_meta` object.
#' @param type What to draw. `"forest"` (default): the studies and the pooled effect of the strongest features;
#'   `"heatmap"`: studies by features; `"volcano"`: pooled estimate against evidence; `"alpha"`: forest plot of the
#'   alpha-diversity metrics; `"beta"`: PERMANOVA R2 of each study; `"influence"`: the pooled effect with each study left out in turn (for features estimated in three or more studies).
#'   `"alpha"` and `"beta"` need `alpha = TRUE` or `beta = TRUE` in the analysis.
#' @param feature The features to draw in a forest, influence or heatmap plot (default: the `top` features with the smallest p-value).
#' @param top Number of features to draw: the strongest ones by p-value (default 3 for forest and influence plots, 25 for the heatmap).
#' @param ... Further arguments for the plot function ([plot_forest()], [plot_heatmap()] or [plot_volcano()]).
#' @rdname meta_analyse
#' @export
plot.mm_meta <- function(x, type = c("forest", "heatmap", "volcano", "alpha", "beta", "influence"), feature = NULL, top = NULL, term = NULL, ...) {
  type <- match.arg(type)
  one_term <- function(tab) {
    terms <- unique(tab$term)
    if (is.null(term)) { if (length(terms) != 1) stop("Several terms were pooled (", paste(terms, collapse = ", "), "); choose one with 'term'.", call. = FALSE); terms } else term
  }
  if (type == "alpha") {
    if (is.null(x$alpha)) stop("Alpha diversity was not analysed; use meta_analyse(..., alpha = TRUE).", call. = FALSE)
    a <- as.data.frame(x$alpha); tm <- one_term(a)
    return(plot_forest(x$alpha_effects, x$alpha, feature = unique(a$feature[a$term == tm]), term = tm, ...))
  }
  if (type == "beta") {
    if (is.null(x$beta_fits)) stop("Beta diversity was not analysed; use meta_analyse(..., beta = TRUE).", call. = FALSE)
    return(.plot_beta(x$beta_fits, x$beta))
  }
  p <- as.data.frame(x$pooled)
  if (type == "volcano") return(plot_volcano(x$pooled, ...))
  term <- one_term(p)
  pt <- p[p$term == term, , drop = FALSE]
  if (type == "heatmap") return(plot_heatmap(x$effects, x$pooled, term = term, features = feature, top = if (is.null(top)) 25 else top, ...))
  if (type == "influence") {       # leaving a study out needs at least two left, so only features estimated in three or more studies qualify
    pt3 <- pt[pt$k >= 3, , drop = FALSE]
    if (!nrow(pt3)) stop("No feature was estimated in three or more studies, so there is nothing to leave out.", call. = FALSE)
    if (is.null(feature)) feature <- utils::head(pt3$feature[order(pt3$p_value)], if (is.null(top)) 3 else top)
    else if (length(bad <- setdiff(feature, pt3$feature))) stop("Leaving a study out needs a feature estimated in at least three studies; not so for ", .show(bad), ".", call. = FALSE)
  }
  if (is.null(feature)) feature <- utils::head(pt$feature[order(pt$p_value)], if (is.null(top)) 3 else top)
  switch(type,
    forest = plot_forest(x$effects, x$pooled, feature = feature, term = term, ...),
    influence = .plot_influence(x, feature, term))
}
