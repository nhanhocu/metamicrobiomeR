# Plots reached through plot(<result of meta_analyse()>, type): diversity, beta diversity, influence of single studies.
# They share the style of the other plots (theme_minimal, size 11, the same colours).

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
