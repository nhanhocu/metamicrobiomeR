#' metamicrobiomeR: differential abundance and meta-analysis of microbiome studies
#'
#' Analyse microbiome studies one at a time and pool the results across studies. See the
#' vignette `getting-started` for the workflow.
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
