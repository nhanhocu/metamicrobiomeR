# Results of other tools as an effect table

#' Effect table from the output of another tool
#'
#' Converts per-study results of other differential-abundance tools into an effect table, so
#' that studies analysed with MaAsLin2, ANCOM-BC2 or LinDA (or any tool, with `tool = "generic"`)
#' can be pooled with [pool_effects()] or combined with [combine_pvalues()]. Every study must be
#' analysed with the same tool, model and transformation: [bind_effects()] refuses to combine
#' different tools. You must state the scale of the estimates (`effect_scale`), because it depends
#' on the tool's settings (for example the MaAsLin2 normalisation and transformation).
#'
#' The expected layouts are those documented for each tool at the time of writing (check them
#' against your version): MaAsLin2's `all_results.tsv` (columns `feature`, `metadata`, `value`,
#' `coef`, `stderr`, `pval`, `N`); ANCOM-BC2's `res` table (columns `taxon` and, for each
#' variable, `lfc_*`, `se_*`, `p_*`); LinDA's table for one variable (row names are the features;
#' columns `log2FoldChange`, `lfcSE`, `pvalue`).
#'
#' @param x The result table (a data.frame, or the path of a tab-separated file).
#' @param tool `"maaslin2"`, `"ancombc2"`, `"linda"` or `"generic"`.
#' @param study Study ID.
#' @param effect_scale Scale of the estimates, e.g. `"log_fold_change"` (required).
#' @param term For `"linda"` the name of the variable (LinDA tables do not contain it). For the
#'   other tools, optionally restricts the result to this term.
#' @param columns For `"generic"`: a named list giving the column names for `feature`, `term`,
#'   `estimate`, `se` and `p_value` (and optionally `n`).
#' @param part,rank Model part and feature rank to record.
#' @return An `mm_effects` object with `engine` set to the tool.
#' @family importers
#' @export
effects_from <- function(x, tool = c("generic", "maaslin2", "ancombc2", "linda"), study, effect_scale, term = NULL, columns = NULL,
                         part = "abundance", rank = NA_character_) {
  tool <- match.arg(tool)
  if (missing(study) || !is.character(study) || length(study) != 1) stop("'study' must be a single study ID.", call. = FALSE)
  if (missing(effect_scale) || !is.character(effect_scale) || length(effect_scale) != 1 || !nzchar(effect_scale)) {
    stop("State the scale of the estimates with 'effect_scale' (for example \"log_fold_change\").", call. = FALSE)
  }
  d <- .read_table(x)
  need <- function(cols) { m <- setdiff(cols, colnames(d)); if (length(m)) stop("The table for tool '", tool, "' lacks the ", .word(length(m), "column"), ": ", paste(m, collapse = ", "), ".", call. = FALSE) }
  out <- switch(tool,
    maaslin2 = {
      need(c("feature", "metadata", "value", "coef", "stderr", "pval"))
      data.frame(feature = d$feature, term = if ("name" %in% colnames(d)) d$name else paste0(d$metadata, d$value), estimate = d$coef,
                 se = d$stderr, p_value = d$pval, n = if ("N" %in% colnames(d)) d$N else NA_integer_, stringsAsFactors = FALSE)
    },
    ancombc2 = {
      need("taxon")
      lf <- grep("^lfc_", colnames(d), value = TRUE); nm <- sub("^lfc_", "", lf); nm <- nm[nm != "(Intercept)"]
      if (!length(nm)) stop("No 'lfc_*' columns found in the ANCOM-BC2 table.", call. = FALSE)
      need(c(paste0("se_", nm), paste0("p_", nm)))
      do.call(rbind, lapply(nm, function(v) data.frame(feature = d$taxon, term = v, estimate = d[[paste0("lfc_", v)]],
                                                       se = d[[paste0("se_", v)]], p_value = d[[paste0("p_", v)]], n = NA_integer_, stringsAsFactors = FALSE)))
    },
    linda = {
      if (is.null(term)) stop("LinDA tables do not name the variable: give it with 'term'.", call. = FALSE)
      need(c("log2FoldChange", "lfcSE", "pvalue"))
      data.frame(feature = rownames(d), term = term, estimate = d$log2FoldChange, se = d$lfcSE, p_value = d$pvalue, n = NA_integer_, stringsAsFactors = FALSE)
    },
    generic = {
      if (is.null(columns) || !all(c("feature", "term", "estimate", "se", "p_value") %in% names(columns))) {
        stop("For tool = \"generic\", give 'columns' with the names of feature, term, estimate, se and p_value.", call. = FALSE)
      }
      need(unlist(columns))
      data.frame(feature = d[[columns$feature]], term = d[[columns$term]], estimate = d[[columns$estimate]], se = d[[columns$se]],
                 p_value = d[[columns$p_value]], n = if (!is.null(columns$n)) d[[columns$n]] else NA_integer_, stringsAsFactors = FALSE)
    })
  if (tool != "linda" && !is.null(term)) out <- out[out$term %in% term, , drop = FALSE]
  if (!nrow(out)) stop("No rows left for that term.", call. = FALSE)
  new_effect_table(study = rep(study, nrow(out)), feature = out$feature, term = out$term, estimate = as.numeric(out$estimate), se = as.numeric(out$se),
                   rank = rank, part = part, p_value = as.numeric(out$p_value), n = as.integer(out$n), engine = tool, effect_scale = effect_scale, converged = TRUE)
}
