# Messages: helpers that keep the wording consistent. See dev/MESSAGES.md for the style guide.

.mm_state <- new.env(parent = emptyenv())
.mm_state$told <- character()

# Informational message for decisions made on the user's behalf and for data that was left out.
# Silenced by options(metamicrobiomeR.quiet = TRUE); `once` shows a message once per session.
.inform <- function(..., once = NULL) {
  if (isTRUE(getOption("metamicrobiomeR.quiet", FALSE))) return(invisible(NULL))
  if (!is.null(once)) {
    if (once %in% .mm_state$told) return(invisible(NULL))
    .mm_state$told <- c(.mm_state$told, once)
  }
  message(...)
}

# The first few IDs, with a count of the rest
.show <- function(x, n = 3) {
  x <- as.character(x)
  if (length(x) <= n) paste(x, collapse = ", ") else paste0(paste(x[seq_len(n)], collapse = ", "), ", ... (", length(x) - n, " more)")
}

.plural <- function(n, one, many = paste0(one, "s")) paste(format(n, big.mark = ","), if (n == 1) one else many)

.pct <- function(x) paste0(format(round(100 * x), trim = TRUE), "%")

# Wrong type of argument: say what was expected, where to get one, and what was received
.check_class <- function(x, cls, arg, how = NULL) {
  if (inherits(x, cls)) return(invisible(x))
  got <- if (is.null(x)) "NULL" else class(x)[1]
  stop("'", arg, "' must be ", if (grepl("^(mm_|[aeiou])", cls)) "an " else "a ", cls, " object",
       if (!is.null(how)) paste0(" (see ", how, ")"), "; got ", if (grepl("^[aeiou]", got)) "an " else "a ", got, ".", call. = FALSE)
}

# Report which IDs differ between two sets, e.g. the sample IDs of the counts and of the sample data
.id_mismatch <- function(label, a_name, a, b_name, b) {
  only_a <- setdiff(a, b); only_b <- setdiff(b, a)
  parts <- c(if (length(only_a)) paste0(length(only_a), " only in ", a_name, " (", .show(only_a), ")"),
             if (length(only_b)) paste0(length(only_b), " only in ", b_name, " (", .show(only_b), ")"))
  if (!length(parts)) return(NULL)
  paste0(label, " do not match between ", a_name, " and ", b_name, ": ", paste(parts, collapse = "; "), ".")
}

# Variables missing from a data frame, with what is available
.missing_vars <- function(miss, avail, where = "the sample data") {
  stop("Variable(s) not in ", where, ": ", paste(miss, collapse = ", "), ". Available: ", .show(avail, 6), ".", call. = FALSE)
}

# Package needed for an optional feature, with how to install it
.need <- function(pkg, what) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  bioc <- pkg %in% c("phyloseq", "SummarizedExperiment", "S4Vectors", "biomformat")
  stop("Package '", pkg, "' is needed for ", what, ". Install it with ",
       if (bioc) paste0("BiocManager::install(\"", pkg, "\")") else paste0("install.packages(\"", pkg, "\")"), ".", call. = FALSE)
}

# A list (or ...) of objects of one class
.check_list_of <- function(xs, cls, how) {
  bad <- which(!vapply(xs, inherits, logical(1), cls))
  if (!length(xs) || length(bad)) {
    got <- if (!length(xs)) "nothing" else paste0("a ", class(xs[[bad[1]]])[1], if (length(xs) > 1) paste0(" (element ", bad[1], ")"))
    stop("Supply ", cls, " objects (see ", how, "); got ", got, ".", call. = FALSE)
  }
  invisible(xs)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
