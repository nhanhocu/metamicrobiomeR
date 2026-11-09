# Checks before pooling: compatible studies, taxon overlap, doubtful estimates

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
#' @family checks before pooling
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
#' @family checks before pooling
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

#' Diagnose the per-study estimates before pooling
#'
#' Pooling weights studies by their precision, so an implausibly small standard error
#' gives one study almost all the weight. Two comparisons are made for each part and term.
#' Between studies, a study's estimate is flagged when its standard error is below `se_ratio`
#' times the typical small standard error of the other studies (their median 5th percentile).
#' Within a study, an estimate is flagged when its standard error is below `within_ratio` times
#' the median standard error of that study, which needs at least `min_estimates` estimates and
#' works with a single study. In fits of real and simulated data, well-identified estimates never
#' fell below 7 per cent of their study's median, while estimates from features with too little
#' information (for example two non-zero samples) fell below 1 per cent.
#' Estimates from fits that did not converge, and (optionally) very large estimates, are flagged too.
#'
#' @param effects An `mm_effects` object.
#' @param se_ratio Flag standard errors below this fraction of the other studies' typical
#'   small standard error. Needs at least two studies.
#' @param within_ratio Flag standard errors below this fraction of the median standard error of
#'   the same study, part and term. `0` turns this comparison off.
#' @param min_estimates Minimum number of estimates in a study, part and term for the
#'   within-study comparison.
#' @param max_abs_estimate Optionally flag estimates with a larger absolute value (only
#'   meaningful on scales with a natural range, such as logits).
#' @param quantile The quantile used as a study's "small standard error" (default 5 per cent).
#' @return An object of class `mm_quality`: `studies` (per study, part and term: the number of
#'   estimates, the proportion converged, the minimum, the chosen quantile and the median of
#'   the standard errors, and how many were flagged) and `rows` (the effect rows with the flag
#'   columns `flag_se` (either comparison), `flag_se_within`, `flag_converged`, `flag_estimate`
#'   and `flagged`).
#' @family checks before pooling
#' @export
check_effects_quality <- function(effects, se_ratio = 0.1, within_ratio = 0.02, min_estimates = 10, max_abs_estimate = NULL, quantile = 0.05) {
  .check_class(effects, "mm_effects", "effects", "fit_study() or bind_effects()")
  check_effects(effects)
  d <- as.data.frame(effects)
  d$flag_se <- FALSE; d$flag_se_within <- .se_within_flags(d, within_ratio, min_estimates); d$flag_converged <- !d$converged
  d$flag_estimate <- if (is.null(max_abs_estimate)) FALSE else !is.na(d$estimate) & abs(d$estimate) > max_abs_estimate
  grp <- paste(d$part, d$term, sep = "\r")
  summ <- list()
  for (g in unique(grp)) {
    ii <- which(grp == g)
    q <- tapply(d$se[ii], d$study[ii], function(s) stats::quantile(s[!is.na(s) & s > 0], quantile, names = FALSE))
    for (s in names(q)) {
      other <- setdiff(names(q), s)
      ref <- if (length(other)) stats::median(q[other], na.rm = TRUE) else NA_real_
      jj <- ii[d$study[ii] == s]
      if (!is.na(ref)) d$flag_se[jj] <- !is.na(d$se[jj]) & d$se[jj] > 0 & d$se[jj] < se_ratio * ref
      se <- d$se[jj]
      summ[[length(summ) + 1]] <- data.frame(study = s, part = d$part[jj][1], term = d$term[jj][1], n = length(jj),
        converged = mean(d$converged[jj]), se_min = suppressWarnings(min(se, na.rm = TRUE)), se_q = unname(q[[s]]),
        se_median = stats::median(se, na.rm = TRUE), reference_se_q = ref, flagged_se = sum(d$flag_se[jj] | d$flag_se_within[jj]),
        flagged_se_within = sum(d$flag_se_within[jj]), flagged_estimate = sum(d$flag_estimate[jj]), stringsAsFactors = FALSE)
    }
  }
  d$flag_se <- d$flag_se | d$flag_se_within
  d$flagged <- d$flag_se | d$flag_converged | d$flag_estimate
  structure(list(studies = do.call(rbind, summ), rows = d,
                 settings = list(se_ratio = se_ratio, within_ratio = within_ratio, min_estimates = min_estimates, max_abs_estimate = max_abs_estimate, quantile = quantile)),
            class = "mm_quality")
}

# Within-study flag: a standard error below `ratio` times the median of the same study, part and term
.se_within_flags <- function(d, ratio = 0.02, min_n = 10) {
  out <- rep(FALSE, nrow(d))
  if (!nrow(d) || is.null(ratio) || ratio <= 0) return(out)
  ok <- !is.na(d$se) & is.finite(d$se) & d$se > 0
  g <- paste(d$study, d$part, d$term, sep = "\r")
  for (k in unique(g[ok])) {
    ii <- which(g == k & ok)
    if (length(ii) < min_n) next
    out[ii] <- d$se[ii] < ratio * stats::median(d$se[ii])
  }
  out
}

#' @export
print.mm_quality <- function(x, ...) {
  s <- x$studies
  cat("<mm_quality> ", nrow(x$rows), " estimates; flagged: ", sum(x$rows$flagged), " (standard error ", sum(x$rows$flag_se),
      ", not converged ", sum(x$rows$flag_converged), ", large estimate ", sum(x$rows$flag_estimate), ")\n", sep = "")
  bad <- s[s$flagged_se > 0 | s$flagged_estimate > 0, , drop = FALSE]
  for (i in seq_len(nrow(bad))) cat(sprintf("  study %s, %s/%s: %d implausibly small standard errors (min %.3g; this study's median %.3g; other studies' typical small value %.3g)\n",
                                            bad$study[i], bad$part[i], bad$term[i], bad$flagged_se[i], bad$se_min[i], bad$se_median[i], bad$reference_se_q[i]))
  invisible(x)
}

#' Drop doubtful estimates before pooling
#'
#' Removes estimates flagged by [check_effects_quality()]: from fits that did not converge
#' and/or with an implausibly small standard error (or very large estimate).
#'
#' @inheritParams check_effects_quality
#' @param drop_nonconverged,drop_implausible_se Which flags to act on.
#' @return The filtered `mm_effects` object, with the removed rows and the reason in the
#'   attribute `dropped`.
#' @family checks before pooling
#' @export
filter_effects <- function(effects, drop_nonconverged = TRUE, drop_implausible_se = TRUE, se_ratio = 0.1, within_ratio = 0.02, max_abs_estimate = NULL) {
  q <- check_effects_quality(effects, se_ratio = se_ratio, within_ratio = within_ratio, max_abs_estimate = max_abs_estimate)
  r <- q$rows
  drop <- (drop_nonconverged & r$flag_converged) | (drop_implausible_se & r$flag_se) | (!is.null(max_abs_estimate) & r$flag_estimate)
  out <- as.data.frame(effects)[!drop, , drop = FALSE]
  class(out) <- class(effects); rownames(out) <- NULL
  dropped <- r[drop, , drop = FALSE]
  dropped$reason <- apply(dropped[, c("flag_converged", "flag_se", "flag_estimate")], 1, function(f) paste(c("not converged", "implausible standard error", "large estimate")[f], collapse = ", "))
  attr(out, "dropped") <- dropped[, c("study", "feature", "part", "term", "estimate", "se", "reason")]
  for (a in c("fit_info", "fit_settings")) if (!is.null(attr(effects, a))) attr(out, a) <- attr(effects, a)
  if (nrow(dropped)) {
    why <- c(if (sum(dropped$flag_converged)) paste(sum(dropped$flag_converged), "not converged"),
             if (sum(dropped$flag_se)) paste(sum(dropped$flag_se), "with an implausibly small standard error"),
             if (sum(dropped$flag_estimate)) paste(sum(dropped$flag_estimate), "with a very large estimate"))
    .inform("Dropped ", nrow(dropped), " of ", nrow(r), " estimates before pooling: ", paste(why, collapse = ", "), " (see attr(, \"dropped\")).")
  }
  out
}
