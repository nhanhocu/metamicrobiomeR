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
#' @export
check_effects <- function(x) {
  if (!inherits(x, "mm_effects")) stop("'x' must be an mm_effects object.", call. = FALSE)
  miss <- setdiff(.mm_effect_cols, colnames(x))
  if (length(miss)) stop("Effect table is missing column(s): ", paste(miss, collapse = ", "), ".", call. = FALSE)
  if (anyNA(x$study) || anyNA(x$feature) || anyNA(x$term) || anyNA(x$part)) {
    stop("'study', 'feature', 'part' and 'term' must not be missing.", call. = FALSE)
  }
  if (!all(x$part %in% .mm_parts)) {
    stop("'part' must be one of: ", paste(.mm_parts, collapse = ", "), ".", call. = FALSE)
  }
  if (anyNA(x$converged)) stop("'converged' must not be missing.", call. = FALSE)
  if (any(x$se < 0, na.rm = TRUE)) stop("'se' must not be negative.", call. = FALSE)
  if (any(x$se == 0, na.rm = TRUE)) warning("Some standard errors are exactly zero.", call. = FALSE)
  if (any(is.infinite(x$estimate)) || any(is.infinite(x$se))) stop("'estimate' and 'se' must be finite or NA.", call. = FALSE)
  key <- paste(x$study, x$feature, x$part, x$term, sep = "\r")
  if (anyDuplicated(key)) stop("Duplicate study/feature/part/term rows in the effect table.", call. = FALSE)
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
#' @return An `mm_effects` object.
#' @export
bind_effects <- function(..., allow_mixed_scale = FALSE) {
  xs <- list(...)
  if (length(xs) == 1 && is.list(xs[[1]]) && !inherits(xs[[1]], "data.frame")) xs <- xs[[1]]
  if (!length(xs) || !all(vapply(xs, inherits, logical(1), "mm_effects"))) {
    stop("Supply mm_effects objects (see new_effect_table()).", call. = FALSE)
  }
  out <- do.call(rbind, lapply(xs, as.data.frame))
  class(out) <- c("mm_effects", "data.frame")
  check_effects(out)
  if (!allow_mixed_scale) {
    grp <- paste(out$part, out$term, sep = " / ")
    nscale <- tapply(out$effect_scale, grp, function(s) length(unique(s[!is.na(s)])))
    bad <- names(nscale)[nscale > 1]
    if (length(bad)) {
      stop("Studies report different effect scales for: ", paste(bad, collapse = "; "),
           ". Pooling needs one scale per term (set allow_mixed_scale = TRUE to override).", call. = FALSE)
    }
  }
  rownames(out) <- NULL
  out
}
