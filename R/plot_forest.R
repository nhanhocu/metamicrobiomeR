# Forest plot of per-study effects and the pooled effect

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
#' @export
plot_forest <- function(effects, pooled = NULL, feature, term = NULL, part = NULL,
                        level = 0.95, show_pi = TRUE, feature_labels = NULL) {
  .check_class(effects, "mm_effects", "effects", "fit_study() or bind_effects()")
  check_effects(effects)
  if (missing(feature) || !length(feature)) stop("'feature' is required.", call. = FALSE)
  miss <- setdiff(feature, effects$feature)
  if (length(miss)) stop("Feature(s) not in the effect table: ", paste(utils::head(miss, 3), collapse = ", "), ".", call. = FALSE)
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
  fmt <- function(x) formatC(x, digits = 2, format = "g")
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
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free_x") +
    ggplot2::labs(x = paste0("Effect (", paste(scale_lab, collapse = ", "), "); ", round(100 * level), "% intervals"), y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0, size = 9.5), axis.text.y = ggplot2::element_text(size = 9.5))
}
