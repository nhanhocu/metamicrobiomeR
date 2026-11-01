# Importers that read files: DADA2, QIIME 2, BIOM, MetaPhlAn and function tables

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

# QIIME 2 metadata file -> data.frame with the sample IDs as row names and typed columns
.read_q2_metadata <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path, ".", call. = FALSE)
  md <- utils::read.delim(path, check.names = FALSE, colClasses = "character", na.strings = c("", "NA"), comment.char = "")
  q2 <- grepl("^#q2:", md[[1]]); types <- NULL
  if (any(q2)) {
    tr <- which(q2 & grepl("^#q2:types", md[[1]]))
    if (length(tr)) types <- stats::setNames(unlist(md[tr[1], -1], use.names = FALSE), colnames(md)[-1])
    md <- md[!q2, , drop = FALSE]
  }
  ids <- md[[1]]; md <- md[, -1, drop = FALSE]
  for (v in colnames(md)) {
    ty <- if (!is.null(types) && v %in% names(types)) types[[v]] else NA_character_
    md[[v]] <- if (identical(ty, "numeric")) {
      x <- suppressWarnings(as.numeric(md[[v]]))
      if (any(is.na(x) & !is.na(md[[v]]))) warning("Column '", v, "' is declared numeric in the metadata but has values that are not numbers.", call. = FALSE)
      x
    } else if (identical(ty, "categorical")) factor(md[[v]]) else utils::type.convert(md[[v]], as.is = TRUE)
  }
  rownames(md) <- ids
  md
}

# Keep the samples that are in both the abundance table and the sample data, in the order of the table, and say what was left out
.align_samples <- function(m, samples, what = "the abundance table") {
  if (!is.data.frame(samples)) return(list(m = m, samples = samples))
  have <- colnames(m); keep <- have[have %in% rownames(samples)]
  if (!length(keep)) {
    stop("None of the samples in ", what, " (", .show(have), ") is in the sample data (", .show(rownames(samples)), "). ",
         "The row names of the sample data must be the sample IDs.", call. = FALSE)
  }
  only_t <- setdiff(have, keep); only_s <- setdiff(rownames(samples), keep)
  if (length(only_s)) .inform(.plural(length(only_s), "sample"), " in the sample data ", if (length(only_s) == 1) "is" else "are", " not in ", what, " and ",
                              if (length(only_s) == 1) "is" else "are", " left out (", .show(only_s), ").")
  if (length(only_t)) .inform(.plural(length(only_t), "sample"), " in ", what, " ", if (length(only_t) == 1) "has" else "have", " no sample data and ",
                              if (length(only_t) == 1) "is" else "are", " left out (", .show(only_t), ").")
  list(m = m[, keep, drop = FALSE], samples = samples[keep, , drop = FALSE])
}

# QIIME 2 artifacts (.qza) are zip files with the data under <uuid>/data/
.qza_member <- function(path, name) {
  if (!file.exists(path)) stop("File not found: ", path, ".", call. = FALSE)
  inside <- utils::unzip(path, list = TRUE)$Name
  hit <- inside[grepl(paste0("(^|/)data/", name, "$"), inside)]
  if (!length(hit)) {
    stop("'", basename(path), "' has no data/", name, ". Is it the right kind of QIIME 2 artifact? It holds: ",
         .show(inside[grepl("/data/", inside)], 4), ".", call. = FALSE)
  }
  dir <- tempfile("qza"); utils::unzip(path, files = hit[1], exdir = dir)
  file.path(dir, hit[1])
}

.read_qza_table <- function(path) {
  file <- .qza_member(path, "feature-table.biom")
  .need("biomformat", "reading a QIIME 2 feature-table artifact (.qza); you can also export it and convert it to text, see ?read_qiime2")
  b <- tryCatch(biomformat::read_biom(file), error = function(e) {
    stop("The feature table in '", basename(path), "' could not be read (", conditionMessage(e), "). QIIME 2 stores tables in the HDF5 format, which needs the ",
         "'rhdf5' package (BiocManager::install(\"rhdf5\")); or export the table and convert it to text, see ?read_qiime2.", call. = FALSE)
  })
  m <- as.matrix(biomformat::biom_data(b)); mode(m) <- "numeric"
  m
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
#'   Samples that are in only one of the abundance table and the sample data are left out, and a message says which.
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
  al <- .align_samples(counts, samples)
  as_study(counts = al$m, taxonomy = tax, samples = al$samples, study = .default_type(study, "asv"))
}

#' Read QIIME 2 output
#'
#' Reads the three things a QIIME 2 project gives you: the feature table, the taxonomy and the
#' sample metadata. The feature table and the taxonomy can be the native artifacts (`.qza`) or the
#' text files you get by exporting them; the metadata is the usual tab-separated QIIME 2 metadata file.
#'
#' To export an artifact instead of reading it directly:
#' ```
#' qiime tools export --input-path table.qza --output-path export
#' biom convert -i export/feature-table.biom -o feature-table.tsv --to-tsv
#' qiime tools export --input-path taxonomy.qza --output-path export    # writes taxonomy.tsv
#' ```
#'
#' The metadata file may have the `#q2:types` row: columns declared `numeric` are read as numbers and
#' columns declared `categorical` as factors (so a batch coded 1, 2, 3 is not treated as a number). Without
#' that row, numeric-looking columns are read as numbers. Empty cells are missing values.
#' Samples that are in only one of the table and the metadata (controls left out of the table, samples
#' without metadata) are left out, and the message says which.
#'
#' @param feature_table The feature table: the path of a `.qza` artifact (needs `biomformat` and, for the
#'   HDF5 format QIIME 2 uses, `rhdf5`), the path of a text file exported with `biom convert --to-tsv`, or a
#'   data.frame (features in rows, samples in columns, the first column holding the feature IDs).
#' @param taxonomy The taxonomy: the path of a `.qza` artifact, the path of the exported `taxonomy.tsv`, or a
#'   data.frame with the columns `Feature ID` and `Taxon`. Optional.
#' @param samples Sample metadata: the path of a QIIME 2 metadata file, or a data.frame with the sample IDs as row names.
#' @param study Named list of study metadata; `feature_type` defaults to `"asv"`.
#' @return An `mm_study` with counts.
#' @examples
#' dir <- system.file("extdata", "qiime2-example", package = "metamicrobiomeR")
#' study <- read_qiime2(file.path(dir, "feature-table.tsv"), file.path(dir, "taxonomy.tsv"),
#'                      file.path(dir, "sample-metadata.tsv"), list(id = "example"))
#' study
#' @family importers
#' @export
read_qiime2 <- function(feature_table, taxonomy = NULL, samples, study) {
  is_path <- function(x) is.character(x) && length(x) == 1
  is_qza <- function(x) is_path(x) && grepl("\\.qza$", x, ignore.case = TRUE)
  if (is_qza(feature_table)) {
    counts <- .read_qza_table(feature_table)
  } else {
    ft <- if (is_path(feature_table) && file.exists(feature_table)) {
      ln <- readLines(feature_table, warn = FALSE)
      ln <- ln[!grepl("^# Constructed from biom file", ln)]
      utils::read.delim(text = ln, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
    } else .read_table(feature_table)
    ids <- as.character(ft[[1]]); ft <- ft[, -1, drop = FALSE]
    counts <- as.matrix(ft); rownames(counts) <- ids; mode(counts) <- "numeric"
  }
  ids <- rownames(counts)
  tax <- NULL
  if (!is.null(taxonomy)) {
    tt <- if (is_qza(taxonomy)) utils::read.delim(.qza_member(taxonomy, "taxonomy.tsv"), check.names = FALSE, stringsAsFactors = FALSE) else .read_table(taxonomy)
    if (!"Taxon" %in% colnames(tt)) stop("'taxonomy' needs a column called 'Taxon' (and the feature IDs in its first column, called 'Feature ID' in QIIME 2 exports).", call. = FALSE)
    tid <- as.character(tt[[if ("Feature ID" %in% colnames(tt)) "Feature ID" else 1]])
    n_no_tax <- sum(!ids %in% tid)
    if (n_no_tax) .inform(n_no_tax, " of ", length(ids), " features have no row in the taxonomy table and are left unassigned.")
    p <- parse_taxonomy(tt$Taxon)[match(ids, tid), , drop = FALSE]
    tax <- p[, .mm_ranks, drop = FALSE]; rownames(tax) <- ids
    tax <- tax[, colSums(!is.na(tax)) > 0, drop = FALSE]
  }
  if (is.character(samples) && length(samples) == 1) samples <- .read_q2_metadata(samples)
  al <- .align_samples(counts, samples)
  as_study(counts = al$m, taxonomy = tax, samples = al$samples, study = .default_type(study, "asv"))
}

#' Read a MetaPhlAn profile
#'
#' Reads a (merged) MetaPhlAn table and keeps the clades of one taxonomic level, whose
#' relative abundances (in percent) sum to the classified fraction of each sample.
#'
#' @param profile Merged MetaPhlAn table (file path or data.frame): a `clade_name` column and
#'   one column per sample. Lines starting with `#` are ignored.
#' @param samples data.frame of sample data; row names are the sample IDs.
#'   Samples that are in only one of the abundance table and the sample data are left out, and a message says which.
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
  al <- .align_samples(vals, samples)
  as_study(relabund = al$m, taxonomy = tax, samples = al$samples, study = .default_type(study, "taxon"),
           relabund_scale = "percent")
}

#' Read a BIOM file
#'
#' @param file Path of a BIOM file (needs the `biomformat` package).
#' @param samples data.frame of sample data with sample IDs as row names. If `NULL`, the
#'   sample metadata stored in the BIOM file is used.
#'   Samples that are in only one of the abundance table and the sample data are left out, and a message says which.
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
  al <- .align_samples(m, samples)
  as_study(counts = al$m, taxonomy = tax, samples = al$samples, study = .default_type(study, "otu"))
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
#'   Samples that are in only one of the abundance table and the sample data are left out, and a message says which.
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
  al <- .align_samples(rel, samples)
  as_study(relabund = al$m, taxonomy = feats, samples = al$samples, study = info)
}
