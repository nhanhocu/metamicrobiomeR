# Importers and exporters for Bioconductor objects (phyloseq, SummarizedExperiment)

#' Convert a phyloseq object to a study
#'
#' @param ps A `phyloseq` object (needs the `phyloseq` package).
#' @param study Named list of study metadata (see [as_study()]).
#' @param abundance `"auto"` (default) treats integer values as counts and other values as
#'   relative abundances; or force `"counts"` or `"relabund"`.
#' @return An `mm_study`. Taxonomy ranks are taken from the taxonomy table (any case).
#' @family importers
#' @export
study_from_phyloseq <- function(ps, study, abundance = c("auto", "counts", "relabund")) {
  .need("phyloseq", "study_from_phyloseq()")
  abundance <- match.arg(abundance)
  m <- methods::as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) m <- t(m)
  tt <- phyloseq::access(ps, "tax_table", errorIfNULL = FALSE)
  tax <- if (is.null(tt)) NULL else .rank_columns(as.data.frame(methods::as(tt, "matrix"), stringsAsFactors = FALSE))
  if (!is.null(tax)) tax <- tax[rownames(m), , drop = FALSE]
  sd <- phyloseq::access(ps, "sam_data", errorIfNULL = FALSE)
  smp <- if (is.null(sd)) data.frame(row.names = colnames(m)) else as.data.frame(methods::as(sd, "data.frame"))
  is_count <- abundance == "counts" || (abundance == "auto" && all(abs(m - round(m)) < 1e-8))
  as_study(counts = if (is_count) m else NULL, relabund = if (!is_count) m else NULL, taxonomy = tax,
           samples = smp, study = .default_type(study, "asv"))
}

#' Convert a SummarizedExperiment or TreeSummarizedExperiment to a study
#'
#' @param se A `SummarizedExperiment` (or `TreeSummarizedExperiment`) with features in rows
#'   (needs the `SummarizedExperiment` package).
#' @param assay Name or number of the assay to use.
#' @param study Named list of study metadata (see [as_study()]).
#' @param abundance As in [study_from_phyloseq()].
#' @return An `mm_study`; the row data supply the taxonomy ranks and the column data the sample data.
#' @family importers
#' @export
study_from_se <- function(se, assay = 1, study, abundance = c("auto", "counts", "relabund")) {
  .need("SummarizedExperiment", "study_from_se()")
  abundance <- match.arg(abundance)
  m <- as.matrix(SummarizedExperiment::assay(se, assay))
  rd <- as.data.frame(SummarizedExperiment::rowData(se)); rownames(rd) <- rownames(m)
  tax <- .rank_columns(rd)
  if (!ncol(tax)) tax <- NULL
  smp <- as.data.frame(SummarizedExperiment::colData(se)); rownames(smp) <- colnames(m)
  is_count <- abundance == "counts" || (abundance == "auto" && all(abs(m - round(m)) < 1e-8))
  as_study(counts = if (is_count) m else NULL, relabund = if (!is_count) m else NULL, taxonomy = tax,
           samples = smp, study = .default_type(study, "asv"))
}

#' Convert a study to a phyloseq or SummarizedExperiment object
#'
#' @param study An `mm_study`.
#' @return A `phyloseq` object (counts if available, else relative abundances), or a
#'   `SummarizedExperiment` with the study's assays.
#' @family importers
#' @name study_export
NULL

#' @rdname study_export
#' @export
as_phyloseq <- function(study) {
  .need("phyloseq", "as_phyloseq()")
  m <- if (!is.null(mm_counts(study))) mm_counts(study) else study$assays$relabund
  f <- mm_features(study)
  ranks <- intersect(.mm_ranks, colnames(f))
  parts <- list(phyloseq::otu_table(m, taxa_are_rows = TRUE), phyloseq::sample_data(mm_samples(study)))
  if (length(ranks)) parts <- c(parts, list(phyloseq::tax_table(as.matrix(f[, ranks, drop = FALSE]))))
  do.call(phyloseq::phyloseq, parts)
}

#' @rdname study_export
#' @export
as_summarized_experiment <- function(study) {
  .need("SummarizedExperiment", "as_summarized_experiment()")
  SummarizedExperiment::SummarizedExperiment(assays = study$assays, rowData = S4Vectors::DataFrame(mm_features(study), check.names = FALSE),
                                             colData = S4Vectors::DataFrame(mm_samples(study), check.names = FALSE))
}
