# Plots: forest, heatmap, volcano, beta diversity, influence

#' Forest plot of per-study and pooled effects
#'
#' Draws, for one or several features, the estimate and confidence interval of
#' every study, the pooled effect as a diamond, and (optionally) the prediction
#' interval. The panel title gives the pooled estimate, its interval and p-value,
#' the heterogeneity (I-squared) and the number of studies. Marker size is
#' proportional to the study's random-effects weight in the pooled result, which is also
#' printed beside each estimate; studies without an estimate for a feature are marked
#' "not estimated".
#'
#' @param effects An `mm_effects` object with the per-study estimates (see [bind_effects()]).
#' @param pooled An `mm_pooled` object from [pool_effects()]. If `NULL`, the effects are
#'   pooled with the default settings of [pool_effects()].
#' @param feature Character vector of one or more feature IDs.
#' @param term,part The model term and part to show. Optional when the selected
#'   features have only one.
#' @param level Confidence level of the per-study intervals (normal approximation).
#' @param show_pi Draw the prediction interval of the pooled effect when it is available
#'   (three or more studies)?
#' @param feature_labels Optional named character vector mapping feature IDs to the
#'   labels shown in the panel titles.
#' @return A `ggplot` object, one panel per feature.
#' @family plots
#' @export
plot_forest <- function(effects, pooled = NULL, feature, term = NULL, part = NULL,
                        level = 0.95, show_pi = TRUE, feature_labels = NULL) {
  .check_class(effects, "mm_effects", "effects", "fit_study() or bind_effects()")
  check_effects(effects)
  if (missing(feature) || !length(feature)) stop("'feature' is required.", call. = FALSE)
  miss <- setdiff(feature, effects$feature)
  if (length(miss)) stop(.cap(.word(length(miss), "feature")), " not found in the effect table: ", .show(miss), ".", call. = FALSE)
  one <- function(x, what) {
    u <- unique(x)
    if (length(u) != 1) stop("Several ", what, "s available (", paste(u, collapse = ", "), "); choose one with '", what, "'.", call. = FALSE)
    u
  }
  sel <- effects[effects$feature %in% feature, , drop = FALSE]
  if (is.null(term)) term <- one(sel$term, "term")
  sel <- sel[sel$term == term, , drop = FALSE]
  if (is.null(part)) part <- one(sel$part, "part")
  sel <- sel[sel$part == part, , drop = FALSE]
  sel <- sel[!is.na(sel$estimate) & !is.na(sel$se), , drop = FALSE]
  if (!nrow(sel)) stop("No estimates for that feature, term and part.", call. = FALSE)
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) stop("'level' must be between 0 and 1.", call. = FALSE)
  zc <- stats::qnorm(0.5 + level / 2)

  if (is.null(pooled)) {
    pooled <- tryCatch(suppressWarnings(pool_effects(effects, term = term, part = part)), error = function(e) NULL)
  } else .check_class(pooled, "mm_pooled", "pooled", "pool_effects()")
  settings <- attr(pooled, "settings")
  pl <- if (is.null(pooled)) NULL else as.data.frame(pooled)[pooled$feature %in% feature & pooled$term == term & pooled$part == part, , drop = FALSE]

  studies <- unique(sel$study); m <- length(studies)
  ypos <- stats::setNames(m:1 + 0.5, studies)
  lab_of <- function(f) if (!is.null(feature_labels) && !is.na(feature_labels[f])) unname(feature_labels[f]) else f
  fmt <- function(x) trimws(formatC(x, digits = 2, format = "g"))
  panel_title <- function(f) {
    r <- pl[pl$feature == f, , drop = FALSE]
    if (is.null(pl) || !nrow(r)) return(paste0(lab_of(f), "\nnot pooled (too few studies)"))
    s <- sprintf("pooled %s [%s, %s]; p = %s; I2 = %s%%; k = %d", fmt(r$estimate), fmt(r$ci_lo), fmt(r$ci_hi),
                 format.pval(r$p_value, digits = 2, eps = 1e-4), formatC(r$I2, digits = 0, format = "f"), r$k)
    if (!is.null(settings)) s <- paste0(s, "\n", settings$method, ", ", settings$test, " test")
    paste0(lab_of(f), "\n", s)
  }
  titles <- vapply(feature, panel_title, character(1))
  panel_of <- function(f) factor(unname(titles[f]), levels = unname(titles))

  d <- data.frame(panel = panel_of(sel$feature), y = unname(ypos[sel$study]), x = sel$estimate,
                  xmin = sel$estimate - zc * sel$se, xmax = sel$estimate + zc * sel$se, stringsAsFactors = FALSE)
  tau2 <- vapply(sel$feature, function(f) if (!is.null(pl) && any(pl$feature == f)) pl$tau2[pl$feature == f][1] else 0, numeric(1))
  wraw <- 1 / (sel$se^2 + tau2)
  share <- stats::ave(wraw, sel$feature, FUN = function(v) v / sum(v))
  d$size <- 1.5 + 3 * stats::ave(share, sel$feature, FUN = function(v) v / max(v))
  d$wlab <- ifelse(share < 0.005, "<1%", sprintf("%.0f%%", 100 * share))
  miss <- do.call(rbind, lapply(feature, function(f) {
    gone <- setdiff(studies, sel$study[sel$feature == f])
    if (!length(gone)) return(NULL)
    data.frame(panel = panel_of(rep(f, length(gone))), y = unname(ypos[gone]))
  }))

  dia <- NULL; pis <- NULL
  if (!is.null(pl) && nrow(pl)) {
    dia <- do.call(rbind, lapply(seq_len(nrow(pl)), function(i) {
      data.frame(panel = panel_of(pl$feature[i]), id = i, x = c(pl$ci_lo[i], pl$estimate[i], pl$ci_hi[i], pl$estimate[i]),
                 y = c(0, 0.25, 0, -0.25))
    }))
    if (show_pi && any(!is.na(pl$pi_lo))) {
      p2 <- pl[!is.na(pl$pi_lo), , drop = FALSE]
      pis <- data.frame(panel = panel_of(p2$feature), x = p2$pi_lo, xend = p2$pi_hi, y = -0.6)
    }
  }

  scale_lab <- unique(sel$effect_scale[!is.na(sel$effect_scale)])
  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_segment(data = d, ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y), colour = "grey25") +
    ggplot2::geom_point(data = d, ggplot2::aes(x = x, y = y, size = size), shape = 15, colour = "grey25") +
    ggplot2::geom_text(data = d, ggplot2::aes(x = Inf, y = y, label = wlab), hjust = 1.05, vjust = -0.7,
                       size = 2.8, colour = "grey40") +
    ggplot2::scale_size_identity()
  if (!is.null(miss)) p <- p + ggplot2::geom_text(data = miss, ggplot2::aes(x = -Inf, y = y, label = "not estimated"),
                                                  hjust = -0.05, size = 3, colour = "grey55", fontface = "italic")
  if (!is.null(dia)) p <- p + ggplot2::geom_polygon(data = dia, ggplot2::aes(x = x, y = y, group = id), fill = "#2C7FB8", colour = "#1B4F72")
  if (!is.null(pis)) p <- p + ggplot2::geom_segment(data = pis, ggplot2::aes(x = x, xend = xend, y = y, yend = y), colour = "#1B4F72", linetype = "dotted", linewidth = 0.7)
  p +
    ggplot2::scale_y_continuous(breaks = c(unname(ypos), 0), labels = c(studies, "Pooled"), limits = c(-0.9, m + 1.1)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.04, 0.2))) +          # room on the right for the weights
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_x") +
    ggplot2::labs(x = paste0("Effect (", paste(.scale_label(scale_lab), collapse = ", "), "); ", round(100 * level), "% intervals"), y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0, size = 9.5), axis.text.y = ggplot2::element_text(size = 9.5))
}

#' Heatmap of per-study and pooled effects
#'
#' One row per feature, one column per study and a final column for the pooled effect. Colour
#' shows the estimate (diverging scale, centred on zero); stars show p-values (`*` below 0.05,
#' `**` below 0.01, `***` below 0.001; for the pooled column the adjusted p-value). Studies
#' without an estimate are grey. The colour range is set by the 95th percentile of the absolute
#' estimates, so a few extreme values do not wash out the rest.
#'
#' @param effects An `mm_effects` object with several studies.
#' @param pooled An `mm_pooled` object; if `NULL`, the effects are pooled with the defaults of
#'   [pool_effects()].
#' @param term,part The model term and part to show (optional when unambiguous).
#' @param features Features to show; by default the `top` features ordered by `order_by`.
#' @param top Number of features when `features` is `NULL`.
#' @param order_by `"p_value"` (smallest pooled p first) or `"estimate"` (largest absolute pooled
#'   estimate first).
#' @param stars Draw the significance stars?
#' @param feature_labels Optional named character vector of labels for the features.
#' @return A `ggplot` object.
#' @family plots
#' @export
plot_heatmap <- function(effects, pooled = NULL, term = NULL, part = NULL, features = NULL, top = 25, order_by = c("p_value", "estimate"),
                         stars = TRUE, feature_labels = NULL) {
  order_by <- match.arg(order_by)
  .check_class(effects, "mm_effects", "effects", "fit_study() or bind_effects()")
  check_effects(effects)
  one <- function(x, what) { u <- unique(x); if (length(u) != 1) stop("Several ", what, "s available (", paste(u, collapse = ", "), "); choose one with '", what, "'.", call. = FALSE); u }
  if (is.null(term)) term <- one(effects$term, "term")
  sel <- effects[effects$term == term, , drop = FALSE]
  if (is.null(part)) part <- one(sel$part, "part")
  sel <- as.data.frame(sel[sel$part == part, , drop = FALSE])
  if (!nrow(sel)) stop("No estimates for that term and part.", call. = FALSE)
  if (is.null(pooled)) pooled <- tryCatch(suppressWarnings(pool_effects(effects, term = term, part = part)), error = function(e) NULL)
  else .check_class(pooled, "mm_pooled", "pooled", "pool_effects()")
  pl <- if (is.null(pooled)) NULL else as.data.frame(pooled)[pooled$term == term & pooled$part == part, , drop = FALSE]
  if (is.null(features)) {
    if (!is.null(pl) && nrow(pl)) {
      pl <- pl[order(if (order_by == "p_value") pl$p_value else -abs(pl$estimate)), ]
      features <- utils::head(pl$feature, top)
    } else {
      mx <- tapply(abs(sel$estimate), sel$feature, max, na.rm = TRUE); features <- names(sort(mx, decreasing = TRUE))[seq_len(min(top, length(mx)))]
    }
  } else if (length(setdiff(features, sel$feature))) stop(.cap(.word(length(setdiff(features, sel$feature)), "feature")), " not found in the effect table: ", .show(setdiff(features, sel$feature)), ".", call. = FALSE)
  studies <- unique(sel$study); cols <- c(studies, "Pooled")
  lab <- function(f) if (!is.null(feature_labels) && !is.na(feature_labels[f])) unname(feature_labels[f]) else f
  grid <- expand.grid(feature = features, column = cols, stringsAsFactors = FALSE)
  m <- match(paste(grid$feature, grid$column), paste(sel$feature, sel$study))
  grid$estimate <- sel$estimate[m]; grid$p <- sel$p_value[m]
  if (!is.null(pl)) {
    mp <- match(grid$feature, pl$feature); isp <- grid$column == "Pooled"
    grid$estimate[isp] <- pl$estimate[mp][isp]; grid$p[isp] <- (if ("p_adjust" %in% colnames(pl)) pl$p_adjust else pl$p_value)[mp][isp]
  }
  grid$label <- ""
  if (stars) grid$label <- ifelse(is.na(grid$p), "", ifelse(grid$p < 0.001, "***", ifelse(grid$p < 0.01, "**", ifelse(grid$p < 0.05, "*", ""))))
  lim <- stats::quantile(abs(grid$estimate), 0.95, na.rm = TRUE); lim <- max(lim, 1e-6)
  lbl <- stats::setNames(make.unique(vapply(features, lab, character(1))), features)        # labels must be unique
  grid$star_col <- ifelse(!is.na(grid$estimate) & abs(grid$estimate) > 0.6 * lim, "white", "grey10")      # readable on the darkest tiles
  grid$x <- factor(grid$column, levels = cols); grid$y <- factor(unname(lbl[grid$feature]), levels = rev(unname(lbl)))
  p <- ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y, fill = estimate)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_vline(xintercept = length(studies) + 0.5, colour = "grey30") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim), na.value = "grey88",
                                  oob = function(x, range) pmin(pmax(x, range[1]), range[2]), name = "Estimate")
  if (stars) p <- p + ggplot2::geom_text(ggplot2::aes(label = label, colour = star_col), size = 3.2, vjust = 0.75) + ggplot2::scale_colour_identity()
  effect_scale <- unique(sel$effect_scale[!is.na(sel$effect_scale)])
  p + ggplot2::labs(x = NULL, y = NULL, title = paste0(term, " (", part, "): ", paste(.scale_label(effect_scale), collapse = ", ")),
                    caption = "Stars: p < 0.05, 0.01, 0.001 (pooled column: adjusted p). Grey: not estimated.") +
    ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(panel.grid = ggplot2::element_blank(), plot.title = ggplot2::element_text(size = 10), plot.caption = ggplot2::element_text(size = 8))
}

#' Volcano plot of pooled results
#'
#' Pooled estimate against the evidence (minus log10 of the p-value), with the features that are
#' significant after adjustment highlighted and the strongest ones labelled. Panels are drawn for
#' each combination of model part and term.
#'
#' @param pooled An `mm_pooled` object.
#' @param alpha Features with an adjusted p-value at or below this are highlighted.
#' @param label_top Number of features to label per panel.
#' @param feature_labels Optional named character vector of labels.
#' @param x_range Range of the horizontal axis. By default it is set from the bulk of the
#'   estimates (the median absolute estimate plus ten times the median absolute deviation, at
#'   least 1), which is robust to many extreme estimates, so that a few of them, typically from
#'   degenerate fits, do not squash the rest; points beyond the range are drawn at the edge as
#'   triangles. Use `Inf` to show everything.
#' @return A `ggplot` object.
#' @family plots
#' @export
plot_volcano <- function(pooled, alpha = 0.05, label_top = 8, feature_labels = NULL, x_range = NULL) {
  .check_class(pooled, "mm_pooled", "pooled", "pool_effects()")
  d <- as.data.frame(pooled)
  d <- d[!is.na(d$estimate) & !is.na(d$p_value), , drop = FALSE]
  if (!nrow(d)) stop("Nothing to plot.", call. = FALSE)
  d$panel <- paste(d$part, d$term, sep = ": "); d$y <- -log10(pmax(d$p_value, 1e-300)); d$x <- d$estimate
  lim <- if (is.null(x_range)) max(1, stats::median(abs(d$estimate)) + 10 * stats::mad(d$estimate)) else x_range
  d$beyond <- is.finite(lim) & abs(d$x) > lim
  if (is.finite(lim)) d$x <- pmax(pmin(d$x, lim), -lim)
  d$significant <- factor(ifelse(!is.na(d$p_adjust) & d$p_adjust <= alpha, "significant", "not significant"), levels = c("not significant", "significant"))
  lab <- function(f) if (!is.null(feature_labels) && !is.na(feature_labels[f])) unname(feature_labels[f]) else f
  d$label <- ""
  for (p in unique(d$panel)) { i <- which(d$panel == p); top <- i[order(d$p_value[i])][seq_len(min(label_top, length(i)))]; d$label[top] <- vapply(d$feature[top], lab, character(1)) }
  ggplot2::ggplot(d, ggplot2::aes(x = x, y = y, colour = significant)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(shape = beyond), alpha = 0.8, size = 2) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.3))) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.8, hjust = -0.1, vjust = -0.4, check_overlap = TRUE, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c("not significant" = "grey55", "significant" = "#B2182B"), drop = FALSE, name = NULL) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::labs(x = "Pooled estimate", y = expression(-log[10](p)),
                  caption = if (any(d$beyond)) paste0(.plural(sum(d$beyond), "estimate"), " beyond the plotted range ", if (sum(d$beyond) == 1) "is" else "are", " drawn at the edge (triangles).")) + ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
}

# Readable names for the effect scales, for axis labels and titles
.scale_label <- function(x) {
  map <- c(logit_mean = "log-odds", asin_sqrt_mean = "arcsine square-root proportion", clr_mean = "centred log-ratio", sd_units = "standard deviations")
  ifelse(x %in% names(map), unname(map[x]), x)
}

.mm_theme <- function() ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), strip.text = ggplot2::element_text(hjust = 0, size = 9.5), axis.text.y = ggplot2::element_text(size = 9.5))

# Beta diversity: PERMANOVA R2 of each study, with the combined result in the panel title
.plot_beta <- function(fits, pooled) {
  d <- do.call(rbind, lapply(fits, as.data.frame))
  pooled <- as.data.frame(pooled)
  d$panel <- vapply(d$term, function(tm) {
    r <- pooled[pooled$term == tm, , drop = FALSE]
    sprintf("%s: combined p = %s, mean R2 = %.2f (k = %d)", tm, format.pval(r$p_value, digits = 2), r$R2, r$k)
  }, character(1))
  d$study <- factor(d$study, levels = rev(unique(d$study)))
  d$label <- paste0("p = ", format.pval(d$p_value, digits = 2))
  mean_line <- data.frame(panel = vapply(pooled$term, function(tm) d$panel[d$term == tm][1], character(1)), x = pooled$R2)
  ggplot2::ggplot(d, ggplot2::aes(x = R2, y = study)) +
    ggplot2::geom_vline(data = mean_line, ggplot2::aes(xintercept = x), colour = "#2C7FB8", linetype = "dotted", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(size = n), shape = 16, colour = "grey25") +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.25, size = 2.8, colour = "grey40") +
    ggplot2::scale_size_continuous(range = c(2, 4.5), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.02, 0.25))) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_x") +
    ggplot2::labs(x = paste0("R2 (share of the variation explained), ", unique(d$distance)[1], " distance"), y = NULL,
                  caption = paste0("Dotted line: mean R2 weighted by sample size. The combined p-value comes from the studies' p-values;\n",
                                   "a permutation p-value cannot be smaller than 1 / (permutations + 1).")) +
    .mm_theme() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), plot.caption = ggplot2::element_text(size = 8))
}

# Influence: the pooled estimate with every study left out in turn
.plot_influence <- function(x, feature, term) {
  p <- as.data.frame(x$pooled)
  rows <- lapply(feature, function(f) {
    all <- p[p$feature == f & p$term == term, , drop = FALSE][1, ]
    l <- leave_one_out(x$effects, feature = f, term = term)
    rbind(data.frame(panel = f, label = "All studies", estimate = all$estimate, ci_lo = all$ci_lo, ci_hi = all$ci_hi, all = TRUE),
          data.frame(panel = f, label = paste("Without", l$omitted), estimate = l$estimate, ci_lo = l$ci_lo, ci_hi = l$ci_hi, all = FALSE))
  })
  d <- do.call(rbind, rows)
  d$label <- factor(d$label, levels = rev(unique(d$label)))
  ref <- d[d$all, , drop = FALSE]
  scale_lab <- unique(x$effects$effect_scale[x$effects$term == term & !is.na(x$effects$effect_scale)])
  ggplot2::ggplot(d, ggplot2::aes(y = label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_vline(data = ref, ggplot2::aes(xintercept = estimate), colour = "#2C7FB8", linetype = "dotted", linewidth = 0.8) +
    ggplot2::geom_segment(ggplot2::aes(x = ci_lo, xend = ci_hi, yend = label, colour = all), linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = estimate, colour = all), size = 2.6) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#2C7FB8", `FALSE` = "grey25"), guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_x") +
    ggplot2::labs(x = paste0("Pooled effect (", paste(.scale_label(scale_lab), collapse = ", "), "), 95% interval"), y = NULL,
                  caption = "Each row pools the studies again with one left out. A row far from the dotted line means that study drives the result.") +
    .mm_theme() + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), plot.caption = ggplot2::element_text(size = 8))
}
