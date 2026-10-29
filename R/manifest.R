# Provenance manifests: what was analysed and how

#' Provenance manifest
#'
#' Summarises what an object contains and how it was produced, so that an analysis can be
#' reported and reproduced: for a study, its metadata and a data summary (dimensions, depth,
#' zeros, ranks, sample variables); for an effect table, the engines, settings and convergence;
#' for a pooled table, the pooling settings. Software versions and the date are included.
#'
#' @param x An `mm_study`, `mm_effects` or `mm_pooled` object.
#' @return A named list of class `mm_manifest`.
#' @export
manifest <- function(x) {
  ver <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  env <- list(date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), R = paste(R.version$major, R.version$minor, sep = "."),
              metamicrobiomeR = ver("metamicrobiomeR"),
              packages = unlist(lapply(c("gamlss", "glmmTMB", "lme4", "metafor", "vegan"), function(p) stats::setNames(ver(p), p))))
  body <- if (inherits(x, "mm_study")) {
    a <- x$assays; ref <- a[[1]]; cnt <- mm_counts(x)
    list(type = "study", study = x$study, features = nrow(ref), samples = ncol(ref), assays = names(a),
         ranks = intersect(.mm_ranks, colnames(x$features)), feature_ranks = if ("rank" %in% colnames(x$features)) sort(unique(x$features$rank)),
         sample_variables = colnames(x$samples),
         depth = if (!is.null(cnt)) stats::setNames(c(min(colSums(cnt)), stats::median(colSums(cnt)), max(colSums(cnt))), c("min", "median", "max")),
         zero_depth_samples = sum(colSums(mm_relabund(x), na.rm = TRUE) == 0),
         proportion_zero = mean((if (!is.null(cnt)) cnt else a[[1]]) == 0))
  } else if (inherits(x, "mm_effects")) {
    fi <- attr(x, "fit_info")
    list(type = "effects", studies = sort(unique(x$study)), engines = sort(unique(x$engine)), effect_scales = sort(unique(x$effect_scale)),
         parts = sort(unique(x$part)), terms = sort(unique(x$term)), rows = nrow(x), features = length(unique(x$feature)),
         converged = if (!is.null(fi)) c(converged = sum(fi$converged), fitted = nrow(fi)) else mean(x$converged),
         fit_settings = attr(x, "fit_settings"))
  } else if (inherits(x, "mm_pooled")) {
    list(type = "pooled", settings = attr(x, "settings"), rows = nrow(x), features = length(unique(x$feature)),
         studies_per_feature = stats::setNames(range(x$k), c("min", "max")),
         features_dominated_by_one_study = sum(x$max_weight > 0.9, na.rm = TRUE))
  } else stop("'x' must be an mm_study, mm_effects or mm_pooled object.", call. = FALSE)
  structure(c(body, list(environment = env)), class = "mm_manifest")
}

#' @export
print.mm_manifest <- function(x, ...) { cat(.manifest_lines(unclass(x)), sep = "\n"); invisible(x) }

.manifest_lines <- function(x, indent = 0) {
  pad <- strrep("  ", indent); out <- character()
  for (n in names(x)) {
    v <- x[[n]]
    if (is.list(v)) out <- c(out, paste0(pad, n, ":"), .manifest_lines(v, indent + 1))
    else out <- c(out, paste0(pad, n, ": ", paste(if (is.null(names(v))) format(v) else paste0(names(v), "=", format(v)), collapse = ", ")))
  }
  out
}

#' Write a manifest to a file
#'
#' @param x An object accepted by [manifest()], or a manifest.
#' @param file Output path.
#' @param format `"text"` (default, readable) or `"json"` (needs the `jsonlite` package).
#' @return The path, invisibly.
#' @export
write_manifest <- function(x, file, format = c("text", "json")) {
  format <- match.arg(format)
  m <- if (inherits(x, "mm_manifest")) x else manifest(x)
  if (format == "json") {
    .need("jsonlite", "write_manifest(format = 'json')")
    writeLines(jsonlite::toJSON(unclass(m), auto_unbox = TRUE, pretty = TRUE, null = "null"), file)
  } else writeLines(.manifest_lines(unclass(m)), file)
  invisible(file)
}
