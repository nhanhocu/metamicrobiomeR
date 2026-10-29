# Effect table: what models return and what pooling consumes (see dev/SPEC.md)

.mm_effect_cols <- c("study", "feature", "rank", "part", "term", "estimate", "se", "df",
                     "p_value", "n", "engine", "effect_scale", "converged")
.mm_parts <- c("abundance", "prevalence", "dispersion", "overall")

#' Create an effect table
#'
#' One row per study x feature x model part x term. Scalars are recycled.
#'
#' @param study,feature,term Character vectors identifying the estimate.
#' @param estimate,se Numeric coefficient and standard error.
#' @param rank Feature rank (`NA` if not applicable).
#' @param part One of "abundance", "prevalence", "dispersion", "overall".
#' @param df Degrees of freedom (`NA` if not defined).
#' @param p_value Within-study p-value (`NA` if not available).
#' @param n Number of samples used.
#' @param engine Name of the engine that produced the estimate.
#' @param effect_scale Scale of `estimate`, e.g. "logit_mean", "log_odds", "smd".
#' @param converged Logical, did the fit converge.
#' @return An object of class `mm_effects` (a data.frame).
#' @family per-study models
#' @export
new_effect_table <- function(study, feature, term, estimate, se,
                             rank = NA_character_, part = "overall", df = NA_real_,
                             p_value = NA_real_, n = NA_integer_, engine = NA_character_,
                             effect_scale = NA_character_, converged = TRUE) {
  x <- data.frame(study = as.character(study), feature = as.character(feature),
                  rank = as.character(rank), part = as.character(part),
                  term = as.character(term), estimate = as.numeric(estimate),
                  se = as.numeric(se), df = as.numeric(df), p_value = as.numeric(p_value),
                  n = as.integer(n), engine = as.character(engine),
                  effect_scale = as.character(effect_scale), converged = as.logical(converged),
                  stringsAsFactors = FALSE)
  class(x) <- c("mm_effects", "data.frame")
  check_effects(x)
}

#' Check an effect table
#'
#' @param x An `mm_effects` object.
#' @return `x`, invisibly, if all checks pass.
#' @family per-study models
#' @export
check_effects <- function(x) {
  .check_class(x, "mm_effects", "x", "fit_study() or bind_effects()")
  miss <- setdiff(.mm_effect_cols, colnames(x))
  if (length(miss)) stop(.cap(.word(length(miss), "column")), " missing from the effect table: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  if (anyNA(x$study) || anyNA(x$feature) || anyNA(x$term) || anyNA(x$part)) {
    stop("'study', 'feature', 'part' and 'term' must not be missing.", call. = FALSE)
  }
  if (!all(x$part %in% .mm_parts)) {
    stop("'part' must be one of: ", paste(.mm_parts, collapse = ", "), ".", call. = FALSE)
  }
  if (anyNA(x$converged)) stop("'converged' must not be missing.", call. = FALSE)
  if (any(x$se < 0, na.rm = TRUE)) stop("'se' must not be negative.", call. = FALSE)
  if (any(x$se == 0, na.rm = TRUE)) {
    n0 <- sum(x$se == 0, na.rm = TRUE)
    warning(.plural(n0, "standard error"), if (n0 == 1) " is" else " are", " exactly zero; such estimates cannot be weighted and are not used when pooling.", call. = FALSE)
  }
  if (any(is.infinite(x$estimate)) || any(is.infinite(x$se))) stop("'estimate' and 'se' must be finite or NA.", call. = FALSE)
  key <- paste(x$study, x$feature, x$part, x$term, sep = "\r")
  if (anyDuplicated(key)) {
    i <- anyDuplicated(key)
    stop("Duplicate rows in the effect table: study '", x$study[i], "', feature '", x$feature[i], "', part '", x$part[i], "', term '", x$term[i],
         "' appears more than once.", call. = FALSE)
  }
  invisible(x)
}

#' Combine effect tables from several studies
#'
#' Stacks per-study effect tables and checks that they can be pooled: no
#' duplicated rows, and one effect scale per model part and term.
#'
#' @param ... `mm_effects` objects, or a single list of them.
#' @param allow_mixed_scale If `FALSE` (default), error when studies report the
#'   same part and term on different effect scales.
#' @param allow_mixed_engine If `FALSE` (default), error when studies report the
#'   same part and term from different engines: estimates of different models are
#'   not comparable even when the scale label is the same.
#' @return An `mm_effects` object.
#' @family per-study models
#' @export
bind_effects <- function(..., allow_mixed_scale = FALSE, allow_mixed_engine = FALSE) {
  xs <- list(...)
  if (length(xs) == 1 && is.list(xs[[1]]) && !inherits(xs[[1]], "data.frame")) xs <- xs[[1]]
  .check_list_of(xs, "mm_effects", "fit_study() or new_effect_table()")
  out <- do.call(rbind, lapply(xs, as.data.frame))
  class(out) <- c("mm_effects", "data.frame")
  check_effects(out)
  .check_comparable(out, allow_mixed_scale, allow_mixed_engine)
  rownames(out) <- NULL
  out
}

# Estimates can only be pooled or combined if they come from one engine and one scale for each part and term
.check_comparable <- function(x, allow_mixed_scale = FALSE, allow_mixed_engine = FALSE) {
  grp <- paste0(x$part, " / ", x$term)
  report <- function(col) {
    u <- tapply(x[[col]], grp, function(s) sort(unique(s[!is.na(s)])))
    bad <- u[lengths(u) > 1]
    if (!length(bad)) return(NULL)
    paste(vapply(names(bad), function(g) paste0(g, " (", paste(bad[[g]], collapse = " vs "), ")"), character(1)), collapse = "; ")
  }
  if (!allow_mixed_scale && !is.null(m <- report("effect_scale"))) {
    stop("Studies report different effect scales for ", m, ". Fit every study with the same engine, or set allow_mixed_scale = TRUE.", call. = FALSE)
  }
  if (!allow_mixed_engine && !is.null(m <- report("engine"))) {
    stop("Studies were fitted with different engines for ", m, ". Fit every study with the same engine, or set allow_mixed_engine = TRUE.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @export
print.mm_effects <- function(x, n = 6, ...) {
  d <- as.data.frame(x)
  if (!all(.mm_effect_cols %in% colnames(d))) { print(d, ...); return(invisible(x)) }      # columns were dropped: plain data frame
  cat("<mm_effects> ", .plural(length(unique(d$study)), "study", "studies"), " (", .show(sort(unique(d$study)), 4), "), ",
      .plural(length(unique(d$feature)), "feature"), ", ", .plural(nrow(d), "row"), "\n", sep = "")
  cat("  engine: ", paste(sort(unique(d$engine)), collapse = ", "), " | scale: ", paste(sort(unique(d$effect_scale)), collapse = ", "),
      " | parts: ", paste(sort(unique(d$part)), collapse = ", "), " | terms: ", .show(sort(unique(d$term)), 4), "\n", sep = "")
  cat("  converged: ", sum(d$converged), " of ", nrow(d), if (any(!d$converged)) "  (see the 'converged' column)", "\n", sep = "")
  if (nrow(d)) {
    cols <- c("study", "feature", "part", "term", "estimate", "se", "p_value")
    cols <- cols[!(cols %in% c("study", "part", "term") & vapply(cols, function(cc) length(unique(d[[cc]])) == 1, logical(1)))]   # constant: in the header
    show <- d[seq_len(min(n, nrow(d))), cols, drop = FALSE]
    print(show, row.names = FALSE, digits = 3)
    if (nrow(d) > n) cat("# ... ", nrow(d) - n, " more rows\n", sep = "")
  }
  invisible(x)
}
