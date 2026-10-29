# Importers: bring data from current tools into a study object

.rank_columns <- function(tax) {
  tax <- as.data.frame(tax, stringsAsFactors = FALSE)
  lc <- tolower(colnames(tax)); lc[lc %in% c("domain", "superkingdom") & !"kingdom" %in% lc] <- "kingdom"
  m <- match(.mm_ranks, lc)
  out <- data.frame(row.names = rownames(tax))
  for (i in which(!is.na(m))) out[[.mm_ranks[i]]] <- .clean_taxon(as.character(tax[[m[i]]]))
  out
}

.read_table <- function(x, ...) {
  if (is.data.frame(x) || is.matrix(x)) return(as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE))
  if (!is.character(x) || length(x) != 1) stop("Supply a file path or a data.frame; got a ", class(x)[1], ".", call. = FALSE)
  if (!file.exists(x)) stop("File not found: ", x, ".", call. = FALSE)
  utils::read.delim(x, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

.default_type <- function(study, type) {
  if (is.null(study$feature_type)) study$feature_type <- type
  study
}

#' Read a DADA2 result
#'
#' @param seqtab Sequence table from `dada2::makeSequenceTable()`: samples in rows,
#'   amplicon sequence variants (ASVs) in columns named by their sequence.
#' @param taxa Taxonomy from `dada2::assignTaxonomy()`: one row per sequence (row names
#'   are the sequences) and the rank columns Kingdom to Species.
#' @param samples data.frame of sample data; row names are the sample IDs.
#' @param study Named list of study metadata (see [as_study()]); `feature_type` defaults to
#'   `"asv"`.
#' @return An `mm_study` with counts. Features are named `ASV1`, `ASV2`, ... and keep their
#'   sequence in the feature data.
#' @family importers
#' @export
read_dada2 <- function(seqtab, taxa = NULL, samples, study) {
  seqtab <- as.matrix(seqtab)
  if (is.null(rownames(seqtab)) || is.null(colnames(seqtab))) stop("'seqtab' needs sample names (rows) and sequences (columns).", call. = FALSE)
  ids <- paste0("ASV", seq_len(ncol(seqtab)))
  tax <- NULL
  if (!is.null(taxa)) {
    taxa <- as.data.frame(taxa, stringsAsFactors = FALSE)
    if (is.null(rownames(taxa)) || !all(colnames(seqtab) %in% rownames(taxa))) {
      if (nrow(taxa) != ncol(seqtab)) stop("'taxa' must have one row per sequence of 'seqtab' (", nrow(taxa), " rows for ", ncol(seqtab), " sequences).", call. = FALSE)
      if (.row_names_info(taxa) > 0) warning("The row names of 'taxa' do not match the sequences of 'seqtab'; assuming the rows are in the same order.", call. = FALSE)
      rownames(taxa) <- colnames(seqtab)
    }
    tax <- .rank_columns(taxa[colnames(seqtab), , drop = FALSE]); rownames(tax) <- ids
  } else tax <- data.frame(row.names = ids)
  tax$sequence <- colnames(seqtab)
  counts <- t(seqtab); rownames(counts) <- ids
  as_study(counts = counts, taxonomy = tax, samples = samples, study = .default_type(study, "asv"))
}

#' Read QIIME 2 exports
#'
#' @param feature_table Feature table exported with `biom convert --to-tsv` (file path or
#'   data.frame): features in rows, samples in columns, the first column holding the feature IDs.
#' @param taxonomy `taxonomy.tsv` exported from a QIIME 2 taxonomy artifact (file path or
#'   data.frame with the columns `Feature ID` and `Taxon`).
#' @param samples Sample metadata: a data.frame with sample IDs as row names, or the path of a
#'   QIIME 2 metadata file (its `#q2:types` row is ignored).
#' @param study Named list of study metadata; `feature_type` defaults to `"asv"`.
#' @return An `mm_study` with counts.
#' @family importers
#' @export
read_qiime2 <- function(feature_table, taxonomy = NULL, samples, study) {
  ft <- if (is.character(feature_table) && length(feature_table) == 1 && file.exists(feature_table)) {
    ln <- readLines(feature_table, warn = FALSE)
    ln <- ln[!grepl("^# Constructed from biom file", ln)]
    utils::read.delim(text = ln, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
  } else .read_table(feature_table)
  ids <- as.character(ft[[1]]); ft <- ft[, -1, drop = FALSE]
  counts <- as.matrix(ft); rownames(counts) <- ids; mode(counts) <- "numeric"
  tax <- NULL
  if (!is.null(taxonomy)) {
    tt <- .read_table(taxonomy)
    if (!all(c("Feature ID", "Taxon") %in% colnames(tt))) stop("'taxonomy' needs the columns 'Feature ID' and 'Taxon'.", call. = FALSE)
    n_no_tax <- sum(!ids %in% tt[["Feature ID"]])
    if (n_no_tax) .inform(n_no_tax, " of ", length(ids), " features have no row in the taxonomy table and are left unassigned.")
    p <- parse_taxonomy(tt$Taxon)[match(ids, tt[["Feature ID"]]), , drop = FALSE]
    tax <- p[, .mm_ranks, drop = FALSE]; rownames(tax) <- ids
    tax <- tax[, colSums(!is.na(tax)) > 0, drop = FALSE]
  }
  if (is.character(samples) && length(samples) == 1) {
    md <- utils::read.delim(samples, check.names = FALSE, stringsAsFactors = FALSE)
    md <- md[!grepl("^#q2:", md[[1]]), , drop = FALSE]
    samples <- md[, -1, drop = FALSE]; rownames(samples) <- md[[1]]
  }
  as_study(counts = counts, taxonomy = tax, samples = samples, study = .default_type(study, "asv"))
}

#' Read a MetaPhlAn profile
#'
#' Reads a (merged) MetaPhlAn table and keeps the clades of one taxonomic level, whose
#' relative abundances (in percent) sum to the classified fraction of each sample.
#'
#' @param profile Merged MetaPhlAn table (file path or data.frame): a `clade_name` column and
#'   one column per sample. Lines starting with `#` are ignored.
#' @param samples data.frame of sample data; row names are the sample IDs.
#' @param study Named list of study metadata; `feature_type` defaults to `"taxon"`.
#' @param level Taxonomic level to keep, default `"species"`.
#' @return An `mm_study` with relative abundances (converted from percent to proportions)
#'   and a `rank` column set to `level`.
#' @family importers
#' @export
read_metaphlan <- function(profile, samples, study, level = "species") {
  level <- match.arg(level, .mm_ranks)
  pt <- if (is.character(profile) && length(profile) == 1 && file.exists(profile)) {
    ln <- readLines(profile, warn = FALSE); ln <- ln[!grepl("^#", ln)]
    utils::read.delim(text = ln, check.names = FALSE, stringsAsFactors = FALSE)
  } else .read_table(profile)
  if (!"clade_name" %in% colnames(pt)) stop("'profile' needs a 'clade_name' column.", call. = FALSE)
  p <- parse_taxonomy(pt$clade_name, sep = "|")
  keep <- !is.na(p$deepest) & p$deepest == level & !grepl("t__", pt$clade_name, fixed = TRUE)
  if (!any(keep)) stop("No clades at level '", level, "' in the profile (deepest levels found: ", .show(unique(stats::na.omit(p$deepest)), 4), ").", call. = FALSE)
  .inform("Kept ", sum(keep), " of ", length(keep), " rows: the clades at level '", level, "'.")
  pt <- pt[keep, , drop = FALSE]; p <- p[keep, , drop = FALSE]
  vals <- as.matrix(pt[, !colnames(pt) %in% c("clade_name", "NCBI_tax_id"), drop = FALSE]); mode(vals) <- "numeric"
  rownames(vals) <- pt$clade_name
  tax <- p[, .mm_ranks, drop = FALSE]; rownames(tax) <- pt$clade_name
  tax <- tax[, colSums(!is.na(tax)) > 0, drop = FALSE]; tax$rank <- level
  as_study(relabund = vals, taxonomy = tax, samples = samples, study = .default_type(study, "taxon"),
           relabund_scale = "percent")
}

#' Read a BIOM file
#'
#' @param file Path of a BIOM file (needs the `biomformat` package).
#' @param samples data.frame of sample data with sample IDs as row names. If `NULL`, the
#'   sample metadata stored in the BIOM file is used.
#' @param study Named list of study metadata (see [as_study()]).
#' @return An `mm_study` with counts and, when the file has taxonomy, rank columns.
#' @family importers
#' @export
read_biom_study <- function(file, samples = NULL, study) {
  .need("biomformat", "read_biom_study()")
  b <- biomformat::read_biom(file)
  m <- as.matrix(biomformat::biom_data(b))
  tax <- NULL
  obs <- biomformat::observation_metadata(b)
  if (!is.null(obs)) {
    lin <- if (is.data.frame(obs)) apply(obs, 1, function(r) paste(r[!is.na(r)], collapse = ";")) else
      vapply(obs, function(o) paste(unlist(o), collapse = ";"), character(1))
    tax <- parse_taxonomy(lin)[, .mm_ranks, drop = FALSE]; rownames(tax) <- rownames(m)
    tax <- tax[, colSums(!is.na(tax)) > 0, drop = FALSE]
  }
  if (is.null(samples)) {
    sm <- biomformat::sample_metadata(b)
    samples <- if (is.null(sm)) data.frame(row.names = colnames(m)) else as.data.frame(sm)
  }
  as_study(counts = m, taxonomy = tax, samples = samples, study = .default_type(study, "otu"))
}

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

#' Read a functional profile (PICRUSt2, HUMAnN, KEGG, ...)
#'
#' Reads a table of functional features (KOs, EC numbers, pathways) per sample and
#' normalises each sample to proportions. Whether the functions were *measured* (shotgun
#' metagenomes, e.g. HUMAnN) or *predicted* from marker genes (e.g. PICRUSt2) must be stated,
#' because the two should not be pooled without flagging the difference.
#'
#' @param table Functional table (file path or data.frame): features in rows, samples in
#'   columns, the first column holding the feature IDs. An optional `description` column is
#'   kept in the feature data.
#' @param samples data.frame of sample data with sample IDs as row names. A column such as
#'   `nsti` (PICRUSt2) can be used with [subset_samples()].
#' @param study Named list of study metadata.
#' @param measured `TRUE` for measured, `FALSE` for predicted functions (required).
#' @param source Name of the tool and version, recorded as the study's `pipeline`.
#' @return An `mm_study` of type `"function"` with relative abundances.
#' @family importers
#' @export
read_function_table <- function(table, samples, study, measured, source = NULL) {
  if (missing(measured) || !is.logical(measured) || length(measured) != 1 || is.na(measured)) {
    stop("State whether the functions were measured (TRUE) or predicted (FALSE) with 'measured'.", call. = FALSE)
  }
  tt <- .read_table(table)
  ids <- as.character(tt[[1]]); tt <- tt[, -1, drop = FALSE]
  desc <- NULL
  if ("description" %in% tolower(colnames(tt))) {
    j <- which(tolower(colnames(tt)) == "description"); desc <- as.character(tt[[j]]); tt <- tt[, -j, drop = FALSE]
  }
  m <- as.matrix(tt); mode(m) <- "numeric"; rownames(m) <- ids
  .inform("Each sample is normalised to proportions, so functions are analysed as relative abundances.", once = "function-normalise")
  tot <- colSums(m)
  rel <- sweep(m, 2, tot, "/"); rel[, tot == 0] <- 0
  feats <- if (is.null(desc)) NULL else data.frame(description = desc, row.names = ids, stringsAsFactors = FALSE)
  info <- study; info$feature_type <- "function"; info$measured <- measured
  if (!is.null(source)) info$pipeline <- source
  as_study(relabund = rel, taxonomy = feats, samples = samples, study = info)
}

#' Subset the samples or features of a study
#'
#' @param study An `mm_study`.
#' @param keep For `subset_samples()`, an expression evaluated in the sample data (for example
#'   `nsti < 0.1`); for `subset_features()`, an expression evaluated in the feature data. A
#'   logical vector is also accepted.
#' @return An `mm_study` with the selected samples or features. Relative abundances that
#'   were stored (not derived from counts) are not renormalised.
#' @family studies
#' @export
subset_samples <- function(study, keep) {
  .check_class(study, "mm_study", "study", "as_study()")
  sel <- eval(substitute(keep), mm_samples(study), parent.frame())
  sel <- .as_selection(sel, ncol(study$assays[[1]]), "samples")
  as_study(counts = .sub(mm_counts(study), j = sel), relabund = .sub(study$assays$relabund, j = sel),
           taxonomy = mm_features(study), samples = mm_samples(study)[sel, , drop = FALSE], study = study$study)
}

#' @rdname subset_samples
#' @export
subset_features <- function(study, keep) {
  .check_class(study, "mm_study", "study", "as_study()")
  sel <- eval(substitute(keep), mm_features(study), parent.frame())
  sel <- .as_selection(sel, nrow(study$assays[[1]]), "features")
  as_study(counts = .sub(mm_counts(study), i = sel), relabund = .sub(study$assays$relabund, i = sel),
           taxonomy = mm_features(study)[sel, , drop = FALSE], samples = mm_samples(study), study = study$study)
}

.as_selection <- function(sel, n, what) {
  if (is.numeric(sel)) sel <- seq_len(n) %in% sel
  if (!is.logical(sel) || length(sel) != n) stop("'keep' must select ", what, " (a logical vector of length ", n, ").", call. = FALSE)
  sel[is.na(sel)] <- FALSE
  if (!any(sel)) stop("No ", what, " left after subsetting.", call. = FALSE)
  sel
}

.sub <- function(m, i = TRUE, j = TRUE) if (is.null(m)) NULL else m[i, j, drop = FALSE]
