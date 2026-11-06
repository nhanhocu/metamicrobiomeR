# Internal helpers for the legacy (1.x) functions.

#' @importFrom dplyr %>% n
#' @importFrom gamlss random
#' @importFrom gamlss.dist BEZI
#' @importFrom stats sd
NULL

utils::globalVariables(c("sequences.per.sample", "pop", "plotvar", "esticat", "pdot", "estimate", "psigcol", "ll", "ul",
                         "esticol", "se", "study", "rel_abund", "taxa", ".", "sumvar", "groupvar", "age.sample", "age.predicted",
                         "studylab"))

.legacy_drop_levels <- function(x) if (is.factor(x)) droplevels(x) else x

# Fits print their summaries and warn once per taxon. The legacy functions report these once, at the end.
.legacy_quiet <- function(expr) {
  noted <- character()
  utils::capture.output(value <- withCallingHandlers(expr, warning = function(w) {
    if (grepl("has not yet converged|vcov has failed|additive.fit convergence", conditionMessage(w))) {
      noted <<- c(noted, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  }))
  attr(value, "legacy_noted") <- unique(noted)
  value
}

.legacy_summary <- function(fitobj) {
  if (inherits(fitobj, "try-error")) return(fitobj)
  .legacy_quiet(try(summary(fitobj, save = TRUE)))
}

.legacy_fit_info <- function(id, testdat, comvar, adjustvar, fitobj) {
  y <- testdat[, id]
  rhs <- paste(setdiff(c(comvar, adjustvar), "random(personid)"), collapse = "+")
  n_coef <- tryCatch(ncol(stats::model.matrix(stats::as.formula(paste("~", rhs)), data = testdat)), error = function(e) NA_integer_)
  data.frame(id = id, n = length(y), n_nonzero = sum(y > 0, na.rm = TRUE), n_coef = n_coef,
             converged = if (inherits(fitobj, "try-error")) NA else isTRUE(fitobj$converged),
             stringsAsFactors = FALSE)
}

.legacy_fit_warnings <- function(info, n.cyc) {
  notconv <- which(!info$converged)
  if (length(notconv)) {
    warning(length(notconv), " of ", nrow(info), " GAMLSS fits did not converge within n.cyc = ", n.cyc,
            " iterations. Estimates and standard errors may change with more iterations. See attr(result, \"fit_info\").", call. = FALSE)
  }
  low <- which(info$n_nonzero <= info$n_coef)
  if (length(low)) {
    warning(length(low), " of ", nrow(info), " taxa have no more non-zero samples than mean-model coefficients, so their estimates and ",
            "standard errors are not reliable: ", paste(utils::head(info$id[low], 5), collapse = ", "),
            if (length(low) > 5) ", ...", ". See attr(result, \"fit_info\").", call. = FALSE)
  }
  invisible(NULL)
}

.legacy_random_warning <- function() {
  warning("GAMLSS standard errors with a subject random effect can be too small, especially for sparse taxa and for variables that ",
          "do not change within a subject. See vignette(\"notes-earlier-results\").", call. = FALSE)
}

# Earlier versions cut adjusted p-value names at the last dot (pval.adjust.sample for age.sample).
# The full name is now used, and the old name stays as an alias so that existing scripts keep working.
.legacy_padjust_alias <- function(tab) {
  if (is.null(tab)) return(tab)
  pcols <- grep("^Pr\\(>\\|t\\|\\)\\.", colnames(tab), value = TRUE)
  full <- sub("^Pr\\(>\\|t\\|\\)", "pval.adjust", pcols)
  old <- sub(".*\\.", "pval.adjust.", pcols)
  for (k in seq_along(pcols)) {
    if (old[k] != full[k] && sum(old == old[k]) == 1 && !old[k] %in% colnames(tab) && full[k] %in% colnames(tab)) tab[, old[k]] <- tab[, full[k]]
  }
  tab
}
