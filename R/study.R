# Study object: the input contract of metamicrobiomeR 2.0 (see dev/SPEC.md)

.mm_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
.mm_feature_types <- c("otu", "asv", "taxon", "function")

.as_assay <- function(x, what) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("'", what, "' must be a numeric matrix or data.frame (features x samples).", call. = FALSE)
  }
  if (is.null(rownames(x)) || is.null(colnames(x))) {
    stop("'", what, "' needs row names (feature IDs) and column names (sample IDs).", call. = FALSE)
  }
  x
}

.align_rows <- function(x, ids, what) {
  if (is.null(rownames(x)) || !setequal(rownames(x), ids) || anyDuplicated(rownames(x))) {
    missing_ids <- setdiff(ids, rownames(x))
    stop("Row names of '", what, "' do not match the assay: ",
         if (length(missing_ids)) paste0("missing ", paste(utils::head(missing_ids, 3), collapse = ", "), " ...")
         else "unexpected or duplicated IDs", call. = FALSE)
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
#'   to the feature IDs) and rank columns (kingdom to species, any case).
#' @param samples data.frame with one row per sample; row names are the sample IDs.
#' @param study Named list of study-level metadata. `id` and `feature_type`
#'   (one of "otu", "asv", "taxon", "function") are required; `region`,
#'   `platform`, `database`, `db_version`, `pipeline`, `measured` and `subject`
#'   (the name of the repeated-measures ID column in `samples`) are recommended.
#' @param relabund_scale Scale of `relabund`: "proportion" (0 to 1) or "percent"
#'   (0 to 100, converted to proportions).
#' @return An object of class `mm_study`.
#' @export
as_study <- function(counts = NULL, relabund = NULL, taxonomy = NULL, samples,
                     study = list(), relabund_scale = c("proportion", "percent")) {
  relabund_scale <- match.arg(relabund_scale)
  if (is.null(counts) && is.null(relabund)) {
    stop("Supply 'counts' and/or 'relabund'.", call. = FALSE)
  }
  if (missing(samples) || !is.data.frame(samples)) {
    stop("'samples' must be a data.frame with one row per sample.", call. = FALSE)
  }
  if (!is.list(study) || is.null(names(study)) && length(study) > 0) {
    stop("'study' must be a named list.", call. = FALSE)
  }
  for (req in c("id", "feature_type")) {
    if (is.null(study[[req]]) || is.na(study[[req]][1]) || !nzchar(study[[req]][1])) {
      stop("study$", req, " is required.", call. = FALSE)
    }
  }
  if (!study$feature_type %in% .mm_feature_types) {
    stop("study$feature_type must be one of: ", paste(.mm_feature_types, collapse = ", "), ".", call. = FALSE)
  }

  assays <- list()
  if (!is.null(counts)) assays$counts <- .as_assay(counts, "counts")
  if (!is.null(relabund)) {
    r <- .as_assay(relabund, "relabund")
    if (relabund_scale == "percent") r <- r / 100
    assays$relabund <- r
  }
  ref <- assays[[1]]
  if (anyDuplicated(rownames(ref)) || anyDuplicated(colnames(ref))) {
    stop("Feature IDs and sample IDs must be unique.", call. = FALSE)
  }
  if (length(assays) == 2) {
    r <- assays$relabund
    if (!setequal(rownames(r), rownames(ref)) || !setequal(colnames(r), colnames(ref))) {
      stop("'counts' and 'relabund' must have the same feature and sample IDs.", call. = FALSE)
    }
    assays$relabund <- r[rownames(ref), colnames(ref), drop = FALSE]
  }

  # samples: one row per assay column, aligned to the assay order
  if (!setequal(rownames(samples), colnames(ref)) || anyDuplicated(rownames(samples))) {
    stop("Row names of 'samples' must match the assay column names (sample IDs).", call. = FALSE)
  }
  samples <- samples[colnames(ref), , drop = FALSE]

  # features: taxonomy in rank columns, aligned to the assay order
  if (is.null(taxonomy)) {
    features <- data.frame(row.names = rownames(ref))
  } else {
    if (!is.data.frame(taxonomy)) stop("'taxonomy' must be a data.frame.", call. = FALSE)
    features <- .align_rows(taxonomy, rownames(ref), "taxonomy")
    rank_cols <- match(.mm_ranks, tolower(colnames(features)))
    names(features)[rank_cols[!is.na(rank_cols)]] <- .mm_ranks[!is.na(rank_cols)]
    present <- intersect(.mm_ranks, colnames(features))
    for (rk in present) {
      v <- as.character(features[[rk]]); v[!nzchar(v)] <- NA_character_
      features[[rk]] <- v
    }
    if (!"rank" %in% colnames(features) && length(present)) {
      deepest <- apply(features[, present, drop = FALSE], 1, function(r) {
        k <- which(!is.na(r)); if (length(k)) present[max(k)] else NA_character_
      })
      features$rank <- unname(deepest)
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
#' @export
check_study <- function(x) {
  if (!inherits(x, "mm_study")) stop("'x' must be an mm_study object.", call. = FALSE)
  a <- x$assays
  ref <- a[[1]]
  if (length(a) == 0) stop("A study needs at least one assay.", call. = FALSE)
  for (nm in names(a)) {
    m <- a[[nm]]
    if (anyNA(m)) stop("Assay '", nm, "' contains missing values.", call. = FALSE)
    if (any(m < 0)) stop("Assay '", nm, "' contains negative values.", call. = FALSE)
    if (!identical(dimnames(m), dimnames(ref))) stop("Assays must have identical dimnames.", call. = FALSE)
  }
  if (!is.null(a$counts) && any(abs(a$counts - round(a$counts)) > 1e-6)) {
    stop("'counts' must be non-negative integers; use 'relabund' for proportions.", call. = FALSE)
  }
  if (!is.null(a$relabund)) {
    grp <- if ("rank" %in% colnames(x$features)) x$features$rank else rep(NA_character_, nrow(ref))
    grp[is.na(grp)] <- "none"
    sums <- vapply(split(seq_len(nrow(ref)), grp),
                   function(i) max(colSums(a$relabund[i, , drop = FALSE])), numeric(1))
    if (any(a$relabund > 1 + 1e-6) || any(sums > 1 + 1e-3)) {
      stop("'relabund' must be proportions in [0, 1] with sums <= 1 within each rank ",
           "(percent data? use relabund_scale = \"percent\").", call. = FALSE)
    }
  }
  if (!identical(rownames(x$features), rownames(ref))) stop("'features' rows must match the assay rows.", call. = FALSE)
  if (!identical(rownames(x$samples), colnames(ref))) stop("'samples' rows must match the assay columns.", call. = FALSE)
  subj <- x$study$subject
  if (!is.na(subj) && !subj %in% colnames(x$samples)) {
    stop("study$subject ('", subj, "') is not a column of 'samples'.", call. = FALSE)
  }
  zero <- if (!is.null(a$counts)) colSums(a$counts) == 0 else colSums(a$relabund) == 0
  if (any(zero)) warning(sum(zero), " sample(s) have zero total abundance.", call. = FALSE)
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
