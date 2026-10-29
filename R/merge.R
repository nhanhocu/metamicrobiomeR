# One-stage (mega-analysis): analyse the studies together, with study as a covariate

#' Merge studies for a one-stage analysis
#'
#' Combines several studies into one study object, so that a single model can be fitted to all
#' samples with the study as a covariate (a "mega-analysis"), as an alternative to analysing
#' each study separately and pooling ([pool_effects()], the "two-stage" approach). The
#' one-stage analysis can be more powerful with small studies, but it assumes the group effect is
#' the same in every study and relies on the study term to absorb batch effects: include `study`
#' in the model formula, and compare with the two-stage result.
#'
#' Features are matched by ID after optional agglomeration to a rank, so use `key = "name"` to
#' match taxa across databases. By default only features present in all studies are kept (then
#' the sequencing depth of the merged study is the depth over the shared features); `"union"` keeps
#' every feature and fills the missing ones with zero, which treats a taxon missing from a
#' study's database as truly absent. If some studies have only proportions, all studies are
#' merged as proportions.
#'
#' @param studies A list of `mm_study` objects (or several as arguments).
#' @param rank Optional rank to agglomerate every study to first (see [agglomerate()]).
#' @param key `"name"` (default) matches taxa by harmonised name, `"lineage"` by lineage.
#' @param features `"intersection"` (default) or `"union"`.
#' @param study_var Name of the sample column that will hold the study (default `"study"`).
#' @return An `mm_study` whose sample IDs are `<study>:<sample>`, with the extra sample column
#'   `study_var`. If every study has a subject column, a combined subject column
#'   `subject_merged` (study plus subject) is set as the study's subject.
#' @export
merge_studies <- function(studies, rank = NULL, key = c("name", "lineage"), features = c("intersection", "union"), study_var = "study") {
  key <- match.arg(key); features <- match.arg(features)
  xs <- .as_study_list(list(studies))
  if (length(xs) < 2) stop("Merging needs at least two studies; got ", length(xs), ".", call. = FALSE)
  ids <- vapply(xs, function(x) x$study$id, character(1))
  if (anyDuplicated(ids)) stop("Study IDs must be unique; '", ids[anyDuplicated(ids)], "' is used twice.", call. = FALSE)
  if (!is.null(rank)) {
    xs <- lapply(xs, function(x) {
      f <- mm_features(x)
      if ("rank" %in% colnames(f) && all(f$rank %in% rank)) x else agglomerate(x, rank, key = key)
    })
  }
  fid <- lapply(xs, function(x) rownames(x$assays[[1]]))
  keep <- if (features == "intersection") Reduce(intersect, fid) else unique(unlist(fid))
  n_all <- length(unique(unlist(fid)))
  if (!length(keep)) stop("No feature is shared by all studies; use features = \"union\", or match the taxa by name with rank = and key = \"name\".", call. = FALSE)
  if (features == "intersection" && length(keep) < 0.5 * n_all) {
    warning("Only ", length(keep), " of ", n_all, " features (", .pct(length(keep) / n_all), ") are present in all studies; consider features = \"union\" or harmonising the taxonomy.", call. = FALSE)
  }
  all_counts <- all(vapply(xs, function(x) !is.null(mm_counts(x)), logical(1)))
  if (!all_counts && any(vapply(xs, function(x) !is.null(mm_counts(x)), logical(1)))) {
    warning("Only some studies have counts: merging relative abundances.", call. = FALSE)
  }
  grab <- function(x) {
    m <- if (all_counts) mm_counts(x) else mm_relabund(x)
    out <- matrix(0, length(keep), ncol(m), dimnames = list(keep, colnames(m)))
    have <- intersect(keep, rownames(m)); out[have, ] <- m[have, , drop = FALSE]
    out
  }
  mats <- lapply(xs, grab)
  for (i in seq_along(mats)) colnames(mats[[i]]) <- paste(ids[i], colnames(mats[[i]]), sep = ":")
  big <- do.call(cbind, mats)
  # sample data: columns shared by all studies, plus the study
  if (study_var %in% unlist(lapply(xs, function(x) colnames(mm_samples(x))))) stop("'", study_var, "' is already a sample column.", call. = FALSE)
  shared <- Reduce(intersect, lapply(xs, function(x) colnames(mm_samples(x))))
  smp <- do.call(rbind, lapply(seq_along(xs), function(i) {
    s <- mm_samples(xs[[i]])[, shared, drop = FALSE]
    for (v in colnames(s)) if (is.factor(s[[v]])) s[[v]] <- as.character(s[[v]])
    rownames(s) <- paste(ids[i], rownames(s), sep = ":"); s[[study_var]] <- ids[i]; s
  }))
  smp[[study_var]] <- factor(smp[[study_var]], levels = ids)
  for (v in setdiff(colnames(smp), study_var)) if (is.character(smp[[v]])) smp[[v]] <- factor(smp[[v]])
  subj <- vapply(xs, function(x) x$study$subject, character(1))
  info <- list(id = paste(ids, collapse = "+"), feature_type = if (length(unique(vapply(xs, function(x) x$study$feature_type, ""))) == 1) xs[[1]]$study$feature_type else "taxon",
               merged_from = ids)
  if (!anyNA(subj)) {
    smp$subject_merged <- unlist(lapply(seq_along(xs), function(i) paste(ids[i], as.character(mm_samples(xs[[i]])[[subj[i]]]), sep = ":")))
    info$subject <- "subject_merged"
  } else if (any(!is.na(subj))) warning("The subject column is missing in some studies; repeated measures are ignored.", call. = FALSE)
  # features: taxonomy from the first study that has each feature
  fcols <- unique(unlist(lapply(xs, function(x) colnames(mm_features(x)))))
  feat <- data.frame(row.names = keep)
  for (cc in fcols) feat[[cc]] <- NA_character_
  for (i in rev(seq_along(xs))) {
    f <- mm_features(xs[[i]]); have <- intersect(keep, rownames(f))
    for (cc in colnames(f)) feat[have, cc] <- as.character(f[have, cc])
  }
  feat <- feat[, colSums(!is.na(feat)) > 0, drop = FALSE]
  out <- as_study(counts = if (all_counts) big else NULL, relabund = if (!all_counts) big else NULL,
                  taxonomy = if (ncol(feat)) feat else NULL, samples = smp, study = info)
  .inform("Merged ", length(xs), " studies: ", .plural(length(keep), "feature"), if (features == "intersection") paste0(" (of ", n_all, ", those present in all studies)")
          else " (missing ones filled with zero)", " in ", .plural(ncol(big), "sample"), ".")
  out
}
