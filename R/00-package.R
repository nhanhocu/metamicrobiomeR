#' metamicrobiomeR: differential abundance and meta-analysis of microbiome studies
#'
#' Analyse microbiome studies one at a time and pool the results across studies. See the
#' vignette `getting-started` for the workflow.
#'
#' @section Where to start:
#' * **One call:** [meta_analyse()] runs the whole workflow (checks, per-study models, quality check,
#'   pooling) and returns an object that prints a summary and plots. Start here.
#' * **A report:** [write_report()] turns the result into a shareable Markdown report with figures.
#' * **Step by step** (full control): build studies with [as_study()] or an importer ([read_dada2()],
#'   [study_from_phyloseq()], ...); fit a model to each with [fit_study()]; combine with [bind_effects()]
#'   and [pool_effects()]; look with [plot_forest()], [plot_heatmap()], [plot_volcano()].
#' * **Before you trust a result:** [check_poolable()], [check_effects_quality()], [pool_sensitivity()],
#'   [leave_one_out()].
#' * **Other tools' results:** [effects_from()] brings MaAsLin2, ANCOM-BC2 or LinDA output into the same
#'   framework.
#'
#' @section Messages:
#' Functions print a short message when they make a decision on your behalf or leave data out
#' (for example samples with missing values, features below the filters, or features pooled from
#' very few studies). Warnings are reserved for things that are likely to affect the results.
#' To silence the messages, use `options(metamicrobiomeR.quiet = TRUE)` or `suppressMessages()`;
#' warnings are not affected. A message that would repeat on every call is shown once per session.
#'
#' @keywords internal
"_PACKAGE"

# Names used in ggplot2 aesthetics (non-standard evaluation)
utils::globalVariables(c("x", "y", "xend", "yend", "xmin", "xmax", "panel", "size", "id", "wlab", "label", "significant", "estimate", "beyond", "star_col", "R2", "ci_lo", "ci_hi", "all", "n", "study", "label", "panel"))
