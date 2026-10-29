# Checks before pooling: are these studies compatible, and how much do their features overlap?

.as_study_list <- function(xs) {
  if (length(xs) == 1 && is.list(xs[[1]]) && !inherits(xs[[1]], "mm_study")) xs <- xs[[1]]
  if (!length(xs) || !all(vapply(xs, inherits, logical(1), "mm_study"))) stop("Supply mm_study objects.", call. = FALSE)
  xs
}

# Names of the taxa a study offers at a rank (harmonised), unassigned excluded
.names_at_rank <- function(study, rank) {
  f <- mm_features(study)
  if (!rank %in% colnames(f)) stop("Study '", study$study$id, "' has no '", rank, "' column in its taxonomy.", call. = FALSE)
  keep <- if ("rank" %in% colnames(f)) f$rank %in% rank else rep(TRUE, nrow(f))
  nm <- harmonise_names(f[[rank]][keep])
  unique(nm[!is.na(nm)])
}

#' Which features are seen in how many studies
#'
#' For an effect table, the studies in which each feature was estimated. For a list of studies
#' and a rank, the taxa each study offers at that rank (names harmonised with
#' [harmonise_names()], unassigned excluded). Features seen in only one study cannot be pooled.
#'
#' @param x An `mm_effects` object, or a list of `mm_study` objects (or several as arguments
#'   via `...` in [check_poolable()]).
#' @param rank For studies: the rank whose taxa are compared.
#' @return A data.frame with `feature`, the number of studies `k` and the study IDs, sorted by
#'   `k`; the logical feature-by-study matrix is in the attribute `membership`.
#' @export
feature_coverage <- function(x, rank = NULL) {
  if (inherits(x, "mm_effects")) {
    d <- unique(as.data.frame(x)[!is.na(x$estimate), c("feature", "study")])
    mem <- table(d$feature, d$study) > 0
  } else {
    xs <- .as_study_list(list(x))
    if (is.null(rank)) stop("Give the 'rank' whose taxa should be compared.", call. = FALSE)
    ids <- vapply(xs, function(s) s$study$id, character(1))
    nm <- lapply(xs, .names_at_rank, rank = rank); all_f <- sort(unique(unlist(nm)))
    mem <- vapply(nm, function(n) all_f %in% n, logical(length(all_f))); dimnames(mem) <- list(all_f, ids)
    if (!length(all_f)) stop("No taxa assigned at rank '", rank, "'.", call. = FALSE)
  }
  out <- data.frame(feature = rownames(mem), k = rowSums(mem),
                    studies = apply(mem, 1, function(r) paste(colnames(mem)[r], collapse = ", ")), stringsAsFactors = FALSE)
  out <- out[order(-out$k, out$feature), ]; rownames(out) <- NULL
  attr(out, "membership") <- mem
  out
}

#' Check that studies can sensibly be pooled
#'
#' Compares the metadata of the studies and reports what is likely to make pooling
#' misleading: mixed feature types, measured and predicted functions in the same analysis,
#' different taxonomy databases, only some studies having counts, very small studies, and,
#' if a rank is given, how many taxa the studies share. Findings are `"warning"` (pooling is
#' likely to mislead without action) or `"note"` (worth knowing; consider a moderator or
#' harmonisation).
#'
#' @param ... `mm_study` objects, or a single list of them.
#' @param rank Optional rank at which to compare the taxa the studies offer (the taxonomy needs
#'   a column for it).
#' @param small_n Studies with fewer samples are flagged.
#' @param strict If `TRUE`, stop when there is a warning.
#' @return An object of class `mm_poolcheck`: `table` (see [study_table()]), `issues` (a
#'   data.frame with `level`, `issue`, `detail`), `coverage` (see [feature_coverage()], or
#'   `NULL`) and `ok` (no warnings).
#' @export
check_poolable <- function(..., rank = NULL, small_n = 10, strict = FALSE) {
  xs <- .as_study_list(list(...))
  st <- study_table(xs)
  iss <- list()
  add <- function(level, issue, detail) iss[[length(iss) + 1]] <<- data.frame(level = level, issue = issue, detail = detail, stringsAsFactors = FALSE)
  ft <- unique(st$feature_type)
  if (length(ft) > 1) add(if ("function" %in% ft) "warning" else "note", "feature types differ", paste(ft, collapse = ", "))
  fn <- st[st$feature_type == "function", , drop = FALSE]
  if (nrow(fn)) {
    if (anyNA(fn$measured)) add("warning", "measured or predicted not stated for some function studies", paste(fn$study[is.na(fn$measured)], collapse = ", "))
    if (length(unique(fn$measured[!is.na(fn$measured)])) > 1) {
      add("warning", "measured and predicted functions are mixed",
          paste0("measured: ", paste(fn$study[fn$measured %in% TRUE], collapse = ", "), "; predicted: ", paste(fn$study[fn$measured %in% FALSE], collapse = ", ")))
    }
  }
  chk <- function(col, label, level = "note") {
    v <- st[[col]]; u <- unique(v[!is.na(v)])
    if (length(u) > 1) add(level, paste(label, "differ"), paste(paste0(st$study, ": ", ifelse(is.na(v), "?", v)), collapse = "; "))
    else if (length(u) == 1 && anyNA(v)) add("note", paste(label, "missing for some studies"), paste(st$study[is.na(v)], collapse = ", "))
  }
  chk("database", "taxonomy databases"); chk("region", "sequencing regions"); chk("platform", "platforms")
  has_counts <- vapply(xs, function(x) !is.null(mm_counts(x)), logical(1))
  if (any(has_counts) && !all(has_counts)) {
    add("note", "only some studies have counts", paste("proportions only:", paste(st$study[!has_counts], collapse = ", "), "(the count-based engine cannot be used for all studies)"))
  }
  small <- st$n_samples < small_n
  if (any(small)) add("warning", paste0("small studies (fewer than ", small_n, " samples)"), paste0(st$study[small], " (n = ", st$n_samples[small], ")", collapse = ", "))
  cov <- NULL
  if (!is.null(rank)) {
    cov <- feature_coverage(xs, rank)
    shared <- sum(cov$k >= 2)
    if (shared == 0) add("warning", paste0("no taxon at rank '", rank, "' is seen in two studies"), "nothing can be pooled; check the taxonomies and harmonise the names")
    else if (shared / nrow(cov) < 0.5) add("note", paste0("few taxa at rank '", rank, "' are shared"), sprintf("%d of %d taxa are seen in two or more studies", shared, nrow(cov)))
  }
  issues <- if (length(iss)) do.call(rbind, iss) else data.frame(level = character(), issue = character(), detail = character())
  out <- structure(list(table = st, issues = issues, coverage = cov, ok = !any(issues$level == "warning")), class = "mm_poolcheck")
  if (strict && !out$ok) stop("The studies are not safe to pool: ", paste(issues$issue[issues$level == "warning"], collapse = "; "), ".", call. = FALSE)
  out
}

#' @export
print.mm_poolcheck <- function(x, ...) {
  cat("<mm_poolcheck> ", nrow(x$table), " studies: ", paste(x$table$study, collapse = ", "), "\n", sep = "")
  if (!nrow(x$issues)) cat("  no issues found\n") else for (i in seq_len(nrow(x$issues))) {
    cat(sprintf("  [%s] %s: %s\n", x$issues$level[i], x$issues$issue[i], x$issues$detail[i]))
  }
  if (!is.null(x$coverage)) cat(sprintf("  taxa: %d, seen in 2+ studies: %d, in all: %d\n", nrow(x$coverage), sum(x$coverage$k >= 2), sum(x$coverage$k == nrow(x$table))))
  invisible(x)
}
