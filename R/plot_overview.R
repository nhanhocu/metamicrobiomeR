# Overview plots of a meta-analysis: heatmap of studies by features, and volcano plot

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
  } else if (length(setdiff(features, sel$feature))) stop("Feature(s) not in the effect table: ", paste(utils::head(setdiff(features, sel$feature), 3), collapse = ", "), ".", call. = FALSE)
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
  grid$x <- factor(grid$column, levels = cols); grid$y <- factor(unname(lbl[grid$feature]), levels = rev(unname(lbl)))
  p <- ggplot2::ggplot(grid, ggplot2::aes(x = x, y = y, fill = estimate)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_vline(xintercept = length(studies) + 0.5, colour = "grey30") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim), na.value = "grey88",
                                  oob = function(x, range) pmin(pmax(x, range[1]), range[2]), name = "Estimate")
  if (stars) p <- p + ggplot2::geom_text(ggplot2::aes(label = label), size = 3.2, vjust = 0.75)
  effect_scale <- unique(sel$effect_scale[!is.na(sel$effect_scale)])
  p + ggplot2::labs(x = NULL, y = NULL, title = paste0(term, " (", part, "): ", paste(effect_scale, collapse = ", ")),
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
                  caption = if (any(d$beyond)) paste0(sum(d$beyond), " estimate(s) beyond the plotted range are drawn at the edge (triangles).")) + ggplot2::theme_minimal(base_size = 11) + ggplot2::theme(legend.position = "bottom")
}
