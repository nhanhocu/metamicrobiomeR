# Study object: the input contract (see dev/SPEC.md)

.as_study_list <- function(xs) {
  if (length(xs) == 1 && is.list(xs[[1]]) && !inherits(xs[[1]], "mm_study")) xs <- xs[[1]]
  .check_list_of(xs, "mm_study", "as_study()")
  xs
}

.mm_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

.mm_feature_types <- c("otu", "asv", "taxon", "function")

.as_assay <- function(x, what) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("'", what, "' must be a numeric matrix or data.frame with features in rows and samples in columns; got ",
         if (is.matrix(x)) "a non-numeric matrix" else paste0("a ", class(x)[1]), ".", call. = FALSE)
  }
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("'", what, "' needs row names (feature IDs) and column names (sample IDs).", call. = FALSE)
  }
  x
}

.align_rows <- function(x, ids, what) {
  if (is.null(rownames(x)) || !setequal(rownames(x), ids) || anyDuplicated(rownames(x))) {
    stop(.id_mismatch("Feature IDs", "the assay", ids, paste0("'", what, "'"), rownames(x)) %||%
           paste0("The row names of '", what, "' must be the feature IDs, each used once."), call. = FALSE)
  }
  x[ids, , drop = FALSE]
}

#' Create a study object
#'
#' A study is one cohort or dataset processed with a single pipeline. It bundles
#' a feature-by-sample abundance matrix, feature annotations (taxonomy in rank
#' columns), sample data and study-level metadata. See `dev/SPEC.md` for the full
#' specification.
#'
#' @param counts Numeric matrix or data.frame of non-negative integer counts,
#'   features in rows and samples in columns, with row and column names.
#' @param relabund Numeric matrix or data.frame of relative abundances with the
#'   same dimnames as `counts` (if both are given). Optional when `counts` is given.
#' @param taxonomy Optional data.frame with one row per feature (row names equal
#'   to the feature IDs) and rank columns (kingdom to species, any case). An
#'   optional `rank` column states the level at which each feature is defined
#'   (a taxonomic rank for agglomerated tables, or e.g. "leaf" for ASVs/OTUs);
#'   it is never inferred from the taxonomy, because a feature with an
#'   unassigned genus can still be a genus-level feature. Relative abundances
#'   must sum to at most 1 within each `rank` group.
#' @param samples data.frame with one row per sample; row names are the sample IDs.
#' @param study Named list of study-level metadata. `id` and `feature_type`
#'   (one of "otu", "asv", "taxon", "function") are required; `region`,
#'   `platform`, `database`, `db_version`, `pipeline`, `measured` and `subject`
#'   (the name of the repeated-measures ID column in `samples`) are recommended.
#' @param relabund_scale Scale of `relabund`: "proportion" (0 to 1) or "percent"
#'   (0 to 100, converted to proportions).
#' @return An object of class `mm_study`.
#' @family studies
#' @export
as_study <- function(counts = NULL, relabund = NULL, taxonomy = NULL, samples,
                     study = list(), relabund_scale = c("proportion", "percent")) {
  relabund_scale <- match.arg(relabund_scale)
  if (is.null(counts) && is.null(relabund)) {
    stop("Supply 'counts' and/or 'relabund'.", call. = FALSE)
  }
  if (missing(samples) || !is.data.frame(samples)) {
    stop("'samples' must be a data.frame with one row per sample and the sample IDs as row names",
         if (!missing(samples)) paste0("; got a ", class(samples)[1]), ".", call. = FALSE)
  }
  if (!is.list(study) || is.null(names(study)) && length(study) > 0) {
    stop("'study' must be a named list.", call. = FALSE)
  }
  study <- study[!vapply(study, is.null, logical(1))]  # NULL entries would delete the defaults
  for (req in c("id", "feature_type")) {
    if (is.null(study[[req]]) || is.na(study[[req]][1]) || !nzchar(study[[req]][1])) {
      stop("study$", req, " is required; for example study = list(id = \"cohort1\", feature_type = \"asv\").", call. = FALSE)
    }
  }
  if (!study$feature_type %in% .mm_feature_types) {
    stop("study$feature_type must be one of: ", paste(.mm_feature_types, collapse = ", "), "; got '", study$feature_type, "'.", call. = FALSE)
  }

  assays <- list()
  if (!is.null(counts)) assays$counts <- .as_assay(counts, "counts")
  if (!is.null(relabund)) {
    r <- .as_assay(relabund, "relabund")
    if (relabund_scale == "percent") r <- r / 100
    assays$relabund <- r
  }
  ref <- assays[[1]]
  if (anyDuplicated(rownames(ref))) stop("Feature IDs must be unique; '", rownames(ref)[anyDuplicated(rownames(ref))], "' is repeated.", call. = FALSE)
  if (anyDuplicated(colnames(ref))) stop("Sample IDs must be unique; '", colnames(ref)[anyDuplicated(colnames(ref))], "' is repeated.", call. = FALSE)
  if (length(assays) == 2) {
    r <- assays$relabund
    msg <- c(.id_mismatch("Feature IDs", "'counts'", rownames(ref), "'relabund'", rownames(r)),
             .id_mismatch("Sample IDs", "'counts'", colnames(ref), "'relabund'", colnames(r)))
    if (length(msg)) stop(paste(msg, collapse = " "), call. = FALSE)
    assays$relabund <- r[rownames(ref), colnames(ref), drop = FALSE]
  }

  # samples: one row per assay column, aligned to the assay order
  if (!setequal(rownames(samples), colnames(ref)) || anyDuplicated(rownames(samples))) {
    stop(.id_mismatch("Sample IDs", "the abundance table", colnames(ref), "'samples'", rownames(samples)) %||%
           "The row names of 'samples' must be the sample IDs, each used once.",
         " The row names of 'samples' must be the sample IDs.", call. = FALSE)
  }
  samples <- samples[colnames(ref), , drop = FALSE]

  # features: taxonomy in rank columns, aligned to the assay order
  if (is.null(taxonomy)) {
    features <- data.frame(row.names = rownames(ref))
  } else {
    .check_class(taxonomy, "data.frame", "taxonomy")
    features <- .align_rows(taxonomy, rownames(ref), "taxonomy")
    lc <- tolower(colnames(features))
    lc[lc %in% c("domain", "superkingdom") & !"kingdom" %in% lc] <- "kingdom"      # SILVA and others call the top rank 'domain'
    rank_cols <- match(.mm_ranks, lc)
    names(features)[rank_cols[!is.na(rank_cols)]] <- .mm_ranks[!is.na(rank_cols)]
    present <- intersect(.mm_ranks, colnames(features))
    other <- setdiff(colnames(features), c("sequence", "description", "rank"))      # columns that are not taxonomy by design
    if (!length(present) && length(other) && !"rank" %in% colnames(features)) {
      warning("'taxonomy' has no rank columns (kingdom, phylum, class, order, family, genus, species; any case, 'domain' counts as kingdom), ",
              "so rank-based functions such as agglomerate() will not work.", call. = FALSE)
    }
    for (rk in present) {
      v <- as.character(features[[rk]]); v[!nzchar(v)] <- NA_character_
      features[[rk]] <- v
    }
  }

  defaults <- list(id = NA_character_, feature_type = NA_character_, region = NA_character_,
                   platform = NA_character_, database = NA_character_, db_version = NA_character_,
                   pipeline = NA_character_, measured = NA, subject = NA_character_)
  x <- structure(list(assays = assays, features = features, samples = samples,
                      study = utils::modifyList(defaults, study)),
                 class = "mm_study")
  check_study(x)
}

#' Check a study object
#'
#' Validates the rules in `dev/SPEC.md`. Errors on violations, warns on soft
#' problems (e.g. samples with zero depth).
#'
#' @param x An `mm_study` object.
#' @return `x`, invisibly, if all checks pass.
#' @family studies
#' @export
check_study <- function(x) {
  .check_class(x, "mm_study", "x", "as_study()")
  a <- x$assays
  if (length(a) == 0) stop("A study needs at least one assay.", call. = FALSE)
  ref <- a[[1]]
  for (nm in names(a)) {
    m <- a[[nm]]
    if (anyNA(m)) stop("Assay '", nm, "' contains ", .plural(sum(is.na(m)), "missing value"), "; replace or remove them.", call. = FALSE)
    if (any(m < 0)) stop("Assay '", nm, "' contains ", .plural(sum(m < 0), "negative value"), ".", call. = FALSE)
    if (!identical(dimnames(m), dimnames(ref))) stop("Assays must have identical dimnames.", call. = FALSE)
  }
  if (!is.null(a$counts) && any(abs(a$counts - round(a$counts)) > 1e-6)) {
    stop("'counts' must be whole numbers (", .plural(sum(abs(a$counts - round(a$counts)) > 1e-6), "value"), " are not); use 'relabund' for proportions.", call. = FALSE)
  }
  if (!is.null(a$relabund)) {
    grp <- if ("rank" %in% colnames(x$features)) x$features$rank else rep(NA_character_, nrow(ref))
    grp[is.na(grp)] <- "none"
    sums <- vapply(split(seq_len(nrow(ref)), grp),
                   function(i) max(colSums(a$relabund[i, , drop = FALSE])), numeric(1))
    if (any(a$relabund > 1 + 1e-6) || any(sums > 1 + 1e-3)) {
      stop("'relabund' must be proportions (0 to 1) that sum to at most 1 within each rank; the largest sum is ",
           format(max(sums), digits = 4), ". Percent data? Use relabund_scale = \"percent\".", call. = FALSE)
    }
  }
  if (!identical(rownames(x$features), rownames(ref))) stop("'features' rows must match the assay rows.", call. = FALSE)
  if (!identical(rownames(x$samples), colnames(ref))) stop("'samples' rows must match the assay columns.", call. = FALSE)
  subj <- x$study$subject
  if (length(subj) == 1 && !is.na(subj) && !subj %in% colnames(x$samples)) {
    stop("study$subject ('", subj, "') is not a column of 'samples'. Available: ", .show(colnames(x$samples), 6), ".", call. = FALSE)
  }
  zero <- if (!is.null(a$counts)) colSums(a$counts) == 0 else colSums(a$relabund) == 0
  if (any(zero)) warning(.plural(sum(zero), "sample"), " ", if (sum(zero) == 1) "has" else "have", " zero total abundance (",
                         .show(colnames(ref)[zero]), "); kept, but left out of models.", call. = FALSE)
  invisible(x)
}

#' Accessors for study objects
#'
#' @param x An `mm_study` object.
#' @return `mm_counts()`: the count matrix (or `NULL`); `mm_relabund()`: relative
#'   abundances (derived from counts if not stored; zero-depth samples are `NA`);
#'   `mm_depth()`: total counts per sample (`NA` if only proportions are stored);
#'   `mm_features()`, `mm_samples()`: the annotation data.frames; `mm_study_info()`:
#'   the study metadata list.
#' @family studies
#' @name mm_accessors
NULL

#' @rdname mm_accessors
#' @export
mm_counts <- function(x) x$assays$counts

#' @rdname mm_accessors
#' @export
mm_relabund <- function(x) {
  if (!is.null(x$assays$relabund)) return(x$assays$relabund)
  d <- colSums(x$assays$counts)
  out <- sweep(x$assays$counts, 2, d, "/")
  out[, d == 0] <- NA_real_
  out
}

#' @rdname mm_accessors
#' @export
mm_depth <- function(x) {
  if (is.null(x$assays$counts)) return(stats::setNames(rep(NA_real_, ncol(x$assays[[1]])), colnames(x$assays[[1]])))
  colSums(x$assays$counts)
}

#' @rdname mm_accessors
#' @export
mm_features <- function(x) x$features

#' @rdname mm_accessors
#' @export
mm_samples <- function(x) x$samples

#' @rdname mm_accessors
#' @export
mm_study_info <- function(x) x$study

#' @export
print.mm_study <- function(x, ...) {
  s <- x$study
  ref <- x$assays[[1]]
  ranks <- intersect(.mm_ranks, colnames(x$features))
  cat("<mm_study> ", s$id, "\n", sep = "")
  cat("  ", nrow(ref), " features x ", ncol(ref), " samples; feature type: ", s$feature_type, "\n", sep = "")
  cat("  assays: ", paste(names(x$assays), collapse = ", "), "\n", sep = "")
  if (length(ranks)) cat("  ranks: ", paste(ranks, collapse = " > "), "\n", sep = "")
  meta <- c(region = s$region, platform = s$platform, database = s$database)
  meta <- meta[!is.na(meta)]
  if (length(meta)) cat("  ", paste(names(meta), meta, sep = ": ", collapse = "; "), "\n", sep = "")
  invisible(x)
}

#' Study-level table for meta-regression
#'
#' Collects the metadata of several studies into one data.frame, to be used as the
#' `study_data` of [meta_regress()] or to check what is being pooled.
#'
#' @param ... `mm_study` objects, or a single list of them.
#' @return A data.frame with one row per study: `study`, `feature_type`, `region`,
#'   `platform`, `database`, `db_version`, `pipeline`, `measured`, `n_samples`
#'   and `n_features`.
#' @family studies
#' @export
study_table <- function(...) {
  xs <- list(...)
  if (length(xs) == 1 && is.list(xs[[1]]) && !inherits(xs[[1]], "mm_study")) xs <- xs[[1]]
  .check_list_of(xs, "mm_study", "as_study()")
  rows <- lapply(xs, function(x) {
    s <- x$study; a <- x$assays[[1]]
    data.frame(study = s$id, feature_type = s$feature_type, region = s$region, platform = s$platform,
               database = s$database, db_version = s$db_version, pipeline = s$pipeline,
               measured = as.logical(s$measured), n_samples = ncol(a), n_features = nrow(a), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (anyDuplicated(out$study)) stop("Study IDs must be unique; '", out$study[anyDuplicated(out$study)], "' is used twice.", call. = FALSE)
  out
}
