# Operations on studies: subset, bring to a taxonomic rank, merge

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

#' Agglomerate the features of a study to a taxonomic rank
#'
#' Sums the abundances of features that share a taxon at the requested rank (for
#' example ASVs into genera), using the rank columns of the feature data. Features
#' without an assignment at that rank are collected into one feature
#' (`unassigned_<rank>`) so that sample totals, and therefore proportions and
#' sequencing depth, stay unchanged; [fit_study()] leaves such features out of the
#' models by default.
#'
#' @param study An `mm_study` with taxonomy in rank columns whose features are not
#'   yet agglomerated (no `rank` column, or `rank = "leaf"`).
#' @param rank One of kingdom, phylum, class, order, family, genus, species.
#' @param key How taxa are identified: `"lineage"` (default) uses the lineage down to the
#'   rank, so equal names in different lineages stay separate; `"name"` uses the (harmonised)
#'   name at the rank only, which is what matching across studies and databases needs.
#' @param unassigned `"keep"` (default) or `"drop"` the unassigned feature. Dropping
#'   changes the sample totals and so the proportions.
#' @param synonyms Passed to [harmonise_names()] when `key = "name"`.
#' @return An `mm_study` with one feature per taxon, a `rank` column, and the rank
#'   columns down to the requested rank.
#' @family studies
#' @export
agglomerate <- function(study, rank, key = c("lineage", "name"), unassigned = c("keep", "drop"), synonyms = NULL) {
  .check_class(study, "mm_study", "study", "as_study()")
  key <- match.arg(key); unassigned <- match.arg(unassigned)
  rank <- match.arg(rank, .mm_ranks)
  feats <- mm_features(study)
  if (!rank %in% colnames(feats)) {
    stop("The feature data has no '", rank, "' column (it has: ", .show(intersect(.mm_ranks, colnames(feats)), 7) %||% "no rank columns",
         "); agglomerate() needs taxonomy ranks, see parse_taxonomy() or the importers.", call. = FALSE)
  }
  if ("rank" %in% colnames(feats) && any(!is.na(feats$rank) & feats$rank != "leaf")) {
    stop("The study is already agglomerated (its features are at rank ", .show(unique(feats$rank), 3), "); pick a rank with fit_study(rank = ...) instead.", call. = FALSE)
  }
  cols <- intersect(.mm_ranks[seq_len(match(rank, .mm_ranks))], colnames(feats))
  tax <- feats[, cols, drop = FALSE]
  nm <- tax[[rank]]
  id <- if (key == "name") harmonise_names(nm, synonyms) else {
    ifelse(is.na(nm), NA_character_, apply(tax, 1, function(r) paste(r[!is.na(r)], collapse = ";")))
  }
  un <- is.na(id)
  if (any(un)) {
    share <- mean(colSums(mm_relabund(study)[un, , drop = FALSE]), na.rm = TRUE)
    if (unassigned == "keep") {
      .inform(sum(un), " of ", length(un), " features (", .pct(share), " of the abundance) have no assignment at '", rank,
              "'; kept as 'unassigned_", rank, "' and left out of models.")
    } else {
      .inform("Dropped ", sum(un), " of ", length(un), " features (", .pct(share), " of the abundance) that are unassigned at '", rank,
              "'; sample totals now cover assigned features only.")
    }
  }
  id[un] <- paste0("unassigned_", rank)
  if (unassigned == "drop") {
    keep <- !un
    if (!any(keep)) stop("No feature is assigned at rank '", rank, "'; check the taxonomy or use unassigned = \"keep\".", call. = FALSE)
  } else keep <- rep(TRUE, length(id))
  agg <- function(m) if (is.null(m)) NULL else rowsum(m[keep, , drop = FALSE], group = id[keep], reorder = FALSE)
  counts <- agg(mm_counts(study)); rel <- agg(study$assays$relabund)
  first <- !duplicated(id[keep])
  newf <- tax[keep, , drop = FALSE][first, , drop = FALSE]
  rownames(newf) <- id[keep][first]
  if (key == "name") newf[, setdiff(cols, rank)] <- NA_character_
  newf[[rank]] <- if (key == "name") rownames(newf) else newf[[rank]]
  newf[grepl("^unassigned_", rownames(newf)), cols] <- NA_character_
  newf$rank <- rank
  info <- study$study; info$feature_type <- "taxon"
  as_study(counts = counts, relabund = rel, taxonomy = newf, samples = mm_samples(study), study = info)
}

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
#' @family studies
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
  levs <- list()                                                  # keep the level order of factors: the first level is the reference level of a model
  for (x in xs) for (v in shared) if (is.factor(mm_samples(x)[[v]])) levs[[v]] <- union(levs[[v]], levels(mm_samples(x)[[v]]))
  smp <- do.call(rbind, lapply(seq_along(xs), function(i) {
    s <- mm_samples(xs[[i]])[, shared, drop = FALSE]
    for (v in colnames(s)) if (is.factor(s[[v]])) s[[v]] <- as.character(s[[v]])
    rownames(s) <- paste(ids[i], rownames(s), sep = ":"); s[[study_var]] <- ids[i]; s
  }))
  smp[[study_var]] <- factor(smp[[study_var]], levels = ids)
  for (v in setdiff(colnames(smp), study_var)) if (is.character(smp[[v]])) smp[[v]] <- if (!is.null(levs[[v]])) factor(smp[[v]], levels = levs[[v]]) else factor(smp[[v]])
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
