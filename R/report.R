#' Write a report of a meta-analysis
#'
#' Writes a Markdown report with the analysis settings, the studies, the checks, the pooled results,
#' the between-study variation, the sensitivity of the result to the pooling settings, diversity
#' results if they were requested, figures, and the software versions. Figures are saved as PNG
#' files in a folder next to the report. It needs no extra packages; Markdown can be opened in any
#' editor or converted to HTML or Word (for example with pandoc or rmarkdown).
#'
#' @param x An `mm_meta` object from [meta_analyse()].
#' @param file Path of the Markdown file.
#' @param top Number of top features to list.
#' @param figures Save a forest plot, a volcano plot and a heatmap?
#' @return The path of the report, invisibly.
#' @family reporting
#' @export
write_report <- function(x, file = "meta-analysis-report.md", top = 10, figures = TRUE) {
  .check_class(x, "mm_meta", "x", "meta_analyse()")
  s <- x$settings; p <- as.data.frame(x$pooled)
  figdir <- paste0(sub("\\.md$", "", file), "_files"); fignames <- character(); skipped <- character()
  md <- function(d) {
    d <- as.data.frame(d)
    for (j in seq_along(d)) if (is.numeric(d[[j]])) d[[j]] <- ifelse(is.na(d[[j]]), "", formatC(d[[j]], digits = 3, format = "g")) else d[[j]] <- ifelse(is.na(d[[j]]), "", as.character(d[[j]]))
    c(paste0("| ", paste(colnames(d), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
      apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  }
  out <- c("# Meta-analysis report", "", paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), " with metamicrobiomeR ", as.character(utils::packageVersion("metamicrobiomeR")), " on R ", as.character(getRversion()), "."), "",
           "## Analysis", "", paste0("- Model: `", s$formula, "`; engine: `", s$engine, "`; rank: ", if (is.null(s$rank)) "as given" else s$rank),
           paste0("- Pooling: ", s$method, " with the ", s$test, " test; estimates flagged as doubtful: ", if (is.null(x$quality)) "not checked" else sum(x$quality$rows$flagged),
                  if (s$quality == "filter") " (dropped before pooling)" else if (s$quality == "warn") " (kept)" else ""), "")
  st <- x$studies; keep <- intersect(c("study", "n_samples", "n_features", "feature_type", "platform", "region", "database"), colnames(st))
  out <- c(out, "## Studies", "", md(st[, keep, drop = FALSE]), "")
  out <- c(out, "## Checks before pooling", "")
  if (is.null(x$checks)) out <- c(out, "Not run.", "") else if (!nrow(x$checks$issues)) out <- c(out, "No issues found.", "") else
    out <- c(out, paste0("- **", x$checks$issues$level, "**: ", x$checks$issues$issue, " (", x$checks$issues$detail, ")"), "")
  top_tab <- p[order(p$p_value), , drop = FALSE]; top_tab <- top_tab[seq_len(min(top, nrow(top_tab))), intersect(c("feature", "term", "k", "estimate", "ci_lo", "ci_hi", "p_value", "p_adjust", "I2"), colnames(top_tab)), drop = FALSE]
  out <- c(out, paste0("## Pooled results (top ", nrow(top_tab), " of ", nrow(p), " features, by p-value)"), "", md(top_tab), "",
           paste0("Adjusted p <= 0.05: ", sum(p$p_adjust <= 0.05, na.rm = TRUE), " of ", nrow(p), " features."), "")
  out <- c(out, "## Between-study variation", "", "```", utils::capture.output(print(summary(x))), "```", "")
  sens <- tryCatch(pool_sensitivity(x$effects, term = s$term, significance = "p_adjust"), error = function(e) NULL)
  if (!is.null(sens)) out <- c(out, "## Does the result depend on the pooling setting?", "", "Features significant at adjusted p <= 0.05 under each setting:", "",
                               md(data.frame(setting = names(sens$n_significant), significant = unname(sens$n_significant))), "",
                               paste0("Significant under all settings: ", length(sens$robust), " of ", nrow(sens$wide), "."), "")
  if (!is.null(x$alpha)) { a <- as.data.frame(x$alpha); out <- c(out, "## Alpha diversity", "", "Pooled standardised difference between groups:", "", md(a[, intersect(c("feature", "k", "estimate", "ci_lo", "ci_hi", "p_value", "I2"), colnames(a))]), "") }
  if (!is.null(x$beta)) out <- c(out, "## Beta diversity", "", "PERMANOVA per study, p-values combined across studies (R2 is a sample-size weighted mean, without an interval):", "", md(x$beta), "")
  if (figures && !requireNamespace("ggplot2", quietly = TRUE)) {
    .inform("ggplot2 is not installed, so the report has no figures (install.packages(\"ggplot2\")).")
    figures <- FALSE
  }
  if (figures) {
    dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
    save <- function(name, type, w, h) tryCatch({ ggplot2::ggsave(file.path(figdir, name), plot(x, type), width = w, height = h, dpi = 110, bg = "white"); fignames <<- c(fignames, name) },
                                                error = function(e) skipped <<- c(skipped, type))
    save("forest.png", "forest", 7, 7); save("volcano.png", "volcano", 6.5, 5); save("heatmap.png", "heatmap", 7, 6)
    out <- c(out, "## Figures", "", if (length(fignames)) paste0("![", sub("\\.png$", "", fignames), "](", basename(figdir), "/", fignames, ")"),
             if (length(skipped)) c("", paste0("Not drawn: ", paste(skipped, collapse = ", "), ". This happens when several terms were pooled; draw one with `plot(x, type, term = ...)`.")), "")
  }
  pk <- c("metafor", "glmmTMB", "gamlss", "lme4", "vegan", "ggplot2")
  ver <- vapply(pk, function(q) if (nzchar(system.file(package = q))) as.character(utils::packageVersion(q)) else "not installed", character(1))     # reads the version without loading the package
  out <- c(out, "## Software", "", paste0("- ", c("R", "metamicrobiomeR", pk), ": ", c(as.character(getRversion()), as.character(utils::packageVersion("metamicrobiomeR")), ver)), "")
  writeLines(out, file)
  invisible(file)
}
