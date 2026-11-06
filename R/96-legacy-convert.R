# Convert 1.x wide taxa tables into study objects.

.legacy_codes <- c(k = "kingdom", p = "phylum", c = "class", o = "order", f = "family", g = "genus", s = "species")

.parse_legacy_names <- function(nms) {
  parts <- strsplit(nms, "\\.(?=[a-z]__)", perl = TRUE)
  tax <- matrix(NA_character_, length(nms), length(.legacy_codes), dimnames = list(nms, unname(.legacy_codes)))
  level <- integer(length(nms)); standard <- logical(length(nms))
  for (i in seq_along(parts)) {
    tk <- parts[[i]]
    last <- tk[length(tk)]
    pad <- regmatches(last, regexpr("(\\.other)+$", last))
    n_pad <- if (length(pad)) lengths(regmatches(pad, gregexpr("\\.other", pad))) else 0L
    tk[length(tk)] <- sub("(\\.other)+$", "", last)
    cd <- substr(tk, 1, 1)
    standard[i] <- all(substr(tk, 2, 3) == "__") && identical(unname(cd), names(.legacy_codes)[seq_along(cd)]) &&
      length(tk) + n_pad <= length(.legacy_codes)
    if (!standard[i]) next
    level[i] <- length(tk) + n_pad
    val <- substring(tk, 4)
    tax[i, seq_along(tk)] <- ifelse(nzchar(val), val, NA_character_)
  }
  list(tax = as.data.frame(tax, stringsAsFactors = FALSE), level = level, standard = standard)
}

#' Convert a 1.x wide table into a study object
#'
#' Moves a table used by the legacy functions, such as [taxa.compare()], into the modern workflow.
#' Samples are in rows, sample variables and taxa are in columns, and taxa have QIIME 1 style names
#' such as `k__bacteria.p__firmicutes`. A trailing `.other` marks an unassigned deeper level.
#'
#' @param tab Wide table with samples in rows and taxa in columns whose names start with `k__`.
#' @param study List of study metadata, as in [as_study()]. It needs `id` and `feature_type`.
#' @param id_col Name of the column holding sample IDs.
#' @return A study object, as returned by [as_study()]. Duplicated sample IDs are made unique. The table
#'   holds proportions only, so the count-based engine cannot be used on it.
#' @seealso [as_study()], [fit_study()], `vignette("migrating-from-1x")`
#' @examples
#' data(taxtab6)
#' s <- legacy_to_study(as.data.frame(taxtab6),
#'                      list(id = "Bangladesh", feature_type = "otu", subject = "personid"))
#' s
#' @export
legacy_to_study <- function(tab, study, id_col = "x.sampleid") {
  tab <- as.data.frame(tab)
  if (!id_col %in% colnames(tab)) stop("Column '", id_col, "' is not in the table. Use 'id_col' to name the sample ID column.", call. = FALSE)
  taxcols <- grep("^k__", colnames(tab))
  if (!length(taxcols)) stop("No taxa columns found. Their names should start with 'k__'.", call. = FALSE)
  p <- .parse_legacy_names(colnames(tab)[taxcols])
  if (!all(p$standard)) stop(sum(!p$standard), " column names do not follow the 1.x pattern, for example '", colnames(tab)[taxcols][which(!p$standard)[1]], "'.", call. = FALSE)
  p$tax$rank <- unname(.legacy_codes)[p$level]
  ids <- make.unique(as.character(tab[[id_col]]))
  m <- t(as.matrix(tab[, taxcols, drop = FALSE])); colnames(m) <- ids
  samp <- tab[, -taxcols, drop = FALSE]; rownames(samp) <- ids
  as_study(relabund = m, taxonomy = p$tax, samples = samp, study = study)
}
