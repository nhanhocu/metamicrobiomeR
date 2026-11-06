# Per-study differential abundance: the engines of fit_study()

.mm_engines <- c("bezi_gamlss", "bezi_glmmtmb", "lm_asin", "lm_clr", "zibb_glmmtmb")

.mm_engine_scale <- c(bezi_gamlss = "logit_mean", bezi_glmmtmb = "logit_mean", lm_asin = "asin_sqrt_mean", lm_clr = "clr_mean",
                      zibb_glmmtmb = "logit_mean")

.zi_se_limit <- 10   # a logit-scale standard error above this means the effect is not identifiable

#' Fit a model to every feature of a study
#'
#' Fits the same model to each feature and returns the estimates as an effect
#' table. Five engines are available.
#'
#' * `"zibb_glmmtmb"`: a zero-inflated beta-binomial model on the counts with the
#'   sequencing depth as the number of trials (`glmmTMB`), a logit link for the mean
#'   abundance and an optional subject random intercept. Needs counts. Returns the
#'   abundance effects (`part = "abundance"`, scale `"logit_mean"`) and, if
#'   `prevalence` is given, effects on the probability of a structural zero
#'   (`part = "prevalence"`, scale `"logit_zero_prob"`; a positive value means the
#'   feature is more often absent).
#' * `"bezi_gamlss"`: the GAMLSS-BEZI model of Ho et al. (2019): a GAMLSS with a zero-inflated
#'   beta (BEZI) family on relative abundance, logit link for the mean. Only the mean
#'   part depends on the covariates. Scale `"logit_mean"`. With a subject random effect
#'   (`random()`), the standard errors can be too small for sparse features and for
#'   covariates that vary between subjects, so `"bezi_glmmtmb"` is preferable for
#'   longitudinal data (see `vignette("notes-earlier-results")`).
#' * `"bezi_glmmtmb"`: the same zero-inflated beta model on relative abundance, fitted
#'   by maximum likelihood with `glmmTMB` (beta family, logit link for the mean, a
#'   logit model for the probability of a zero) and an optional subject random
#'   intercept, whose variance is estimated by Laplace approximation rather than by
#'   the penalised fitting that GAMLSS uses. Does not need counts. Scale
#'   `"logit_mean"`, so its estimates are on the same scale as `"bezi_gamlss"`.
#'   `prevalence` gives effects on the probability of a zero, as for `"zibb_glmmtmb"`.
#' * `"lm_asin"`: a linear model on the arcsine square root of the relative abundance
#'   (the transformation used by MaAsLin), a linear mixed model with a subject random
#'   intercept when `longitudinal` (p-values from the normal approximation),
#'   otherwise `lm()` with t-based p-values. Scale `"asin_sqrt_mean"`.
#' * `"lm_clr"`: the same linear (mixed) model on the centred log-ratio transform of the
#'   features at the analysed rank, as a compositional sensitivity analysis. A pseudocount
#'   replaces zeros. Scale `"clr_mean"`.
#'
#' Estimates from different engines are on different scales and are not comparable;
#' [bind_effects()] refuses to combine them.
#'
#' Features are kept if they are present (abundance above zero) in at least
#' `min_prevalence` of all samples and their mean relative abundance exceeds
#' `min_mean`. Samples with a missing value in the model are dropped for that feature.
#'
#' A feature that is non-zero in only a few samples gives an estimate that is
#' essentially the difference between those few values, often with a tiny standard
#' error that has no inferential meaning. Such features are therefore not fitted:
#' a feature needs more non-zero samples than the mean model has coefficients
#' (see `min_nonzero`). The features left out are listed in the attribute
#' `low_support`.
#'
#' @param study An `mm_study` object.
#' @param formula One-sided formula with the covariates, e.g. `~ group + age`.
#'   Variables are looked up in the sample data.
#' @param engine Model engine (see Details), default `"bezi_gamlss"`.
#' @param rank Optional character vector of feature ranks to analyse (matched
#'   against the `rank` column of the feature data). `NULL` analyses all features.
#' @param longitudinal Add a subject-level random intercept? Defaults to `TRUE`
#'   when `study$subject` is set. The subject is the sample-data column named in
#'   `study$subject`.
#' @param min_prevalence,min_mean Feature filters (see Details).
#' @param min_nonzero Minimum number of samples, among those used in the model, in
#'   which a feature must be non-zero to be fitted. `NULL` (default) uses the number
#'   of coefficients of the mean model (intercept included) plus one. Use `0` to fit
#'   every feature that passes the other filters.
#' @param n_cycles (`bezi_gamlss` only) Maximum number of outer iterations of the GAMLSS
#'   fitting algorithm (`n.cyc` in `gamlss::gamlss.control()`). The default is 100.
#'   The GAMLSS default is 20, at which many fits on real data do not converge.
#'   Fits that reach the limit are flagged `converged = FALSE` and a warning reports how many.
#' @param exclude_unassigned Drop features that are unassigned at their own rank
#'   (for example a genus-level feature with no genus, such as
#'   `...f__lachnospiraceae.other`)? Such features lump together unrelated
#'   organisms. Default `TRUE`; use `FALSE` to keep them.
#' @param prevalence (`zibb_glmmtmb` and `bezi_glmmtmb`) One-sided formula for the covariates of the
#'   probability of a structural zero, e.g. `~ group`. `NULL` (default) fits an
#'   intercept-only zero-inflation part, so only abundance effects are reported.
#'   Prevalence effects can be poorly identified; unidentifiable ones are flagged
#'   `converged = FALSE`.
#' @param zero_inflation (`zibb_glmmtmb` and `bezi_glmmtmb`) Include a zero-inflation part? It is
#'   dropped automatically for features without zeros. `FALSE` fits a plain
#'   beta-binomial model (`zibb_glmmtmb`), or a beta model in which zeros are replaced by
#'   half the smallest non-zero value of the feature (`bezi_glmmtmb`).
#' @param pseudocount (`lm_clr` only) Value added to counts (default 0.5) or, for
#'   proportions, to relative abundances (default half the smallest positive value).
#' @param cores Number of processor cores for fitting the features in parallel (default 1, or
#'   `options(metamicrobiomeR.cores = )`). Uses forking, so it is available on Linux and macOS but not on
#'   Windows, where one core is used. Worthwhile for the count model and GAMLSS, which take about a second
#'   per feature; the linear engines are fast anyway. Results do not depend on `cores`.
#' @return An `mm_effects` object with one row per feature, model part and covariate
#'   term. The attribute `fit_info` is a data.frame with one row per fitted feature
#'   (samples used, non-zero samples, convergence, number of warnings, and any error
#'   message), and the attribute `low_support` lists the features not fitted because
#'   they had too few non-zero samples.
#' @family per-study models
#' @export
fit_study <- function(study, formula, engine = "bezi_gamlss", rank = NULL,
                      longitudinal = NULL, min_prevalence = 0.05, min_mean = 5e-05,
                      n_cycles = 100, exclude_unassigned = TRUE, prevalence = NULL,
                      zero_inflation = TRUE, pseudocount = NULL, cores = getOption("metamicrobiomeR.cores", 1L),
                      min_nonzero = NULL) {
  .check_class(study, "mm_study", "study", "as_study()")
  engine <- match.arg(engine, .mm_engines)
  if (!inherits(formula, "formula") || length(formula) != 2) {
    stop("'formula' must be one-sided, e.g. ~ group + age.", call. = FALSE)
  }
  if (!is.null(prevalence) && (!inherits(prevalence, "formula") || length(prevalence) != 2)) {
    stop("'prevalence' must be a one-sided formula, e.g. ~ group.", call. = FALSE)
  }
  if (!is.null(min_nonzero) && (!is.numeric(min_nonzero) || length(min_nonzero) != 1 || is.na(min_nonzero) || min_nonzero < 0)) {
    stop("'min_nonzero' must be NULL or a single number of at least 0; got ", format(min_nonzero)[1], ".", call. = FALSE)
  }
  if (!is.numeric(cores) || length(cores) != 1 || is.na(cores) || cores < 1) stop("'cores' must be a single number of at least 1; got ", format(cores)[1], ".", call. = FALSE)
  cores <- as.integer(cores)
  if (cores > 1 && .Platform$OS.type == "windows") {
    .inform("Parallel fitting is not available on Windows; using 1 core.", once = "windows-cores")
    cores <- 1L
  }
  info <- mm_study_info(study)
  sdat <- mm_samples(study)
  vars <- unique(c(all.vars(formula), if (!is.null(prevalence)) all.vars(prevalence)))
  miss <- setdiff(vars, colnames(sdat))
  if (length(miss)) .missing_vars(miss, colnames(sdat))
  subject <- info$subject
  defaulted <- is.null(longitudinal)
  if (defaulted) longitudinal <- !is.na(subject)
  if (defaulted && longitudinal) {
    .inform("Using '", subject, "' as a subject random intercept (set longitudinal = FALSE to ignore it).", once = "longitudinal-default")
  }
  if (longitudinal && engine %in% c("lm_asin", "lm_clr")) .need("lme4", "a subject random effect with the linear engines (longitudinal = FALSE avoids it)")
  if (longitudinal && is.na(subject)) {
    stop("A longitudinal model needs study$subject, the sample column that identifies the subject (see as_study()).", call. = FALSE)
  }
  switch(engine, bezi_gamlss = .need("gamlss", "the 'bezi_gamlss' engine"),
         zibb_glmmtmb = .need("glmmTMB", "the 'zibb_glmmtmb' engine"),
         bezi_glmmtmb = .need("glmmTMB", "the 'bezi_glmmtmb' engine"))
  if (engine == "bezi_gamlss" && longitudinal) {
    .inform("GAMLSS standard errors with a subject random effect can be too small for sparse features and for covariates that vary ",
            "between subjects; engine = \"bezi_glmmtmb\" fits the same model with a mixed-model likelihood.", once = "gamlss-random")
  }
  counts <- mm_counts(study)
  if (engine == "zibb_glmmtmb" && is.null(counts)) {
    stop("The 'zibb_glmmtmb' engine needs counts; this study only has relative abundances.", call. = FALSE)
  }

  d0 <- sdat[, vars, drop = FALSE]
  for (v in vars) if (is.factor(d0[[v]])) d0[[v]] <- droplevels(d0[[v]])
  if (longitudinal) d0$.subject <- factor(sdat[[subject]])

  ra <- mm_relabund(study)
  feats <- mm_features(study)
  idx <- seq_len(nrow(ra))
  if (!is.null(rank)) {
    if (!"rank" %in% colnames(feats)) stop("The feature data has no 'rank' column, so 'rank' cannot select features; agglomerate() creates one.", call. = FALSE)
    idx <- which(feats$rank %in% rank)
    if (!length(idx)) stop("No features with rank ", paste0("'", rank, "'", collapse = ", "), "; the study has: ", .show(sort(unique(feats$rank)), 6), ".", call. = FALSE)
  }
  n_rank <- length(idx); n_unassigned <- 0L
  if (exclude_unassigned && "rank" %in% colnames(feats)) {
    unassigned <- vapply(idx, function(i) {
      r <- feats$rank[i]
      !is.na(r) && r %in% colnames(feats) && is.na(feats[[r]][i])
    }, logical(1))
    n_unassigned <- sum(unassigned)
    idx <- idx[!unassigned]
  }
  composition <- idx                                   # the features that make up the composition (for CLR)
  n_all <- ncol(ra)
  keep <- vapply(idx, function(i) {
    x <- ra[i, ]
    sum(!is.na(x) & x > 0) >= min_prevalence * n_all && mean(x, na.rm = TRUE) > min_mean
  }, logical(1))
  n_filtered <- sum(!keep)
  idx <- idx[keep]
  if (!length(idx)) {
    warning("No feature passed the filters (present in at least ", .pct(min_prevalence), " of the samples and mean relative abundance above ",
            format(min_mean), "); lower 'min_prevalence' or 'min_mean'.", call. = FALSE)
  } else if (n_unassigned + n_filtered > 0) {
    why <- c(if (n_unassigned) paste(n_unassigned, "unassigned"), if (n_filtered) paste(n_filtered, "below the prevalence and abundance filters"))
    .inform("Fitting ", length(idx), " of ", n_rank, " features (", paste(why, collapse = ", "), ").")
  }
  n_cc <- sum(stats::complete.cases(d0))
  if (n_cc < nrow(d0)) .inform(nrow(d0) - n_cc, " of ", nrow(d0), " samples have missing values in the model variables and are left out.")

  # Low support: a feature needs more non-zero samples (among those used) than the mean model has coefficients
  terms_f <- attr(stats::terms(formula), "term.labels")
  n_coef <- tryCatch(ncol(stats::model.matrix(stats::reformulate(c("1", terms_f)), stats::na.omit(d0[, vars, drop = FALSE]))), error = function(e) NA_integer_)
  need_nz <- if (is.null(min_nonzero)) (if (is.na(n_coef)) 0 else n_coef + 1) else min_nonzero
  low_support <- data.frame(feature = character(), n_nonzero = integer(), stringsAsFactors = FALSE)
  if (length(idx) && need_nz > 0) {
    used <- stats::complete.cases(d0)
    if (engine == "zibb_glmmtmb" && !is.null(mm_depth(study))) used <- used & mm_depth(study) > 0
    nz <- vapply(idx, function(i) sum(!is.na(ra[i, used]) & ra[i, used] > 0), numeric(1))
    low <- nz < need_nz
    if (any(low)) {
      low_support <- data.frame(feature = rownames(ra)[idx[low]], n_nonzero = as.integer(nz[low]), stringsAsFactors = FALSE)
      warning(.plural(sum(low), "feature"), " not fitted: non-zero in fewer than ", need_nz, " of the samples used (the mean model has ",
              if (is.na(n_coef)) "an unknown number of" else n_coef, " coefficients), so ", if (sum(low) == 1) "its estimate" else "their estimates",
              " would rest on a few values; see attr(, \"low_support\"), or set 'min_nonzero'.", call. = FALSE)
      idx <- idx[!low]
    }
  }

  Y <- switch(engine,
    bezi_gamlss = ra,
    bezi_glmmtmb = ra,
    zibb_glmmtmb = ra,
    lm_asin = asin(sqrt(ra)),
    lm_clr = .clr_matrix(study, composition, pseudocount))
  depth <- mm_depth(study)

  re_term <- if (engine == "bezi_gamlss") "random(.subject)" else "(1 | .subject)"
  f <- stats::reformulate(c(terms_f, if (longitudinal) re_term), response = ".y")
  # random() must be found without gamlss being attached: evaluate the formula in its namespace
  if (engine == "bezi_gamlss") environment(f) <- asNamespace("gamlss")

  batch <- NULL          # linear engines: one factorisation for all features when they share a design (identical results, far faster)
  if (engine %in% c("lm_asin", "lm_clr") && !longitudinal && length(idx) && n_cc == nrow(d0) && isTRUE(getOption("metamicrobiomeR.batch_lm", TRUE))) {
    batch <- .lm_batch(stats::reformulate(terms_f), d0, Y[rownames(ra)[idx], , drop = FALSE])
  }
  fit_one <- function(k) {
    i <- idx[k]; fid <- rownames(ra)[i]
    d <- d0
    d$.y <- Y[fid, ]
    d$.nz <- ra[fid, ] > 0                                  # non-zero on the original scale, whatever the engine's transform
    if (engine == "zibb_glmmtmb") {
      d$.y <- counts[fid, ]; d$.fail <- depth - d$.y
      d <- d[depth > 0, , drop = FALSE]
    }
    if (engine %in% c("bezi_gamlss", "bezi_glmmtmb")) d$.y[which(d$.y >= 1)] <- 0.9999
    d <- stats::na.omit(d)
    res <- switch(engine,
      bezi_gamlss = .fit_bezi_gamlss(f, d, n_cycles),
      lm_asin = .fit_lm(f, d, longitudinal, .mm_engine_scale[["lm_asin"]]),
      lm_clr = .fit_lm(f, d, longitudinal, .mm_engine_scale[["lm_clr"]]),
      zibb_glmmtmb = .fit_zibb(terms_f, d, longitudinal, prevalence, zero_inflation),
      bezi_glmmtmb = .fit_bezi_glmmtmb(terms_f, d, longitudinal, prevalence, zero_inflation))
    list(log = data.frame(feature = fid, n = nrow(d), n_nonzero = sum(d$.nz), converged = res$converged, n_warnings = length(res$warnings),
                          first_warning = if (length(res$warnings)) res$warnings[1] else NA_character_,
                          singular = isTRUE(res$singular), error = res$error, stringsAsFactors = FALSE),
         row = if (!is.null(res$rows) && nrow(res$rows)) {
           rk <- if ("rank" %in% colnames(feats)) feats$rank[i] else NA_character_
           cbind(feature = fid, rank = rk, res$rows, n = nrow(d), stringsAsFactors = FALSE)
         })
  }
  todo <- if (is.null(batch)) seq_along(idx) else which(!batch$ok)        # features that need lm() themselves (or another engine)
  fits <- if (cores > 1 && length(todo) > 1) parallel::mclapply(todo, fit_one, mc.cores = min(cores, length(todo))) else lapply(todo, fit_one)
  for (j in seq_along(fits)) if (inherits(fits[[j]], "try-error")) {          # a worker died: report it as a failed fit
    fits[[j]] <- list(log = data.frame(feature = rownames(ra)[idx[todo[j]]], n = NA_integer_, n_nonzero = NA_integer_, converged = FALSE, n_warnings = 0L, first_warning = NA_character_,
                                       singular = FALSE, error = trimws(as.character(fits[[j]])), stringsAsFactors = FALSE), row = NULL)
  }
  rows <- lapply(fits, `[[`, "row"); log <- lapply(fits, `[[`, "log")
  tab <- do.call(rbind, rows)
  fit_info <- do.call(rbind, log)
  if (!is.null(batch) && any(batch$ok)) {          # the batch features: one table for all of them, then back into feature order
    fids <- rownames(ra)[idx]; ok_k <- which(batch$ok); terms <- setdiff(rownames(batch$B), "(Intercept)"); nt <- length(terms); nb <- length(ok_k)
    rk_all <- if ("rank" %in% colnames(feats)) feats$rank[idx] else rep(NA_character_, length(idx))
    if (nt) {
      tab_b <- data.frame(feature = rep(fids[ok_k], each = nt), rank = rep(rk_all[ok_k], each = nt), term = rep(terms, times = nb), part = "abundance",
                          estimate = as.vector(batch$B[terms, ok_k, drop = FALSE]), se = as.vector(batch$se[terms, ok_k, drop = FALSE]),
                          p_value = as.vector(batch$pv[terms, ok_k, drop = FALSE]), df = batch$rdf, effect_scale = .mm_engine_scale[[engine]],
                          converged = TRUE, n = nrow(d0), stringsAsFactors = FALSE)
      tab <- rbind(tab, tab_b)
      tab <- tab[order(match(tab$feature, fids)), , drop = FALSE]
    }
    nz_b <- vapply(idx[ok_k], function(i) sum(ra[i, ] > 0, na.rm = TRUE), numeric(1))
    fit_info <- rbind(fit_info, data.frame(feature = fids[ok_k], n = nrow(d0), n_nonzero = nz_b, converged = TRUE, n_warnings = 0L, first_warning = NA_character_,
                                           singular = FALSE, error = NA_character_, stringsAsFactors = FALSE))
    fit_info <- fit_info[order(match(fit_info$feature, fids)), , drop = FALSE]; rownames(fit_info) <- NULL
  }
  if (is.null(tab)) tab <- data.frame(feature = character(), rank = character(), term = character(),
                                      part = character(), estimate = numeric(), se = numeric(), p_value = numeric(),
                                      df = numeric(), effect_scale = character(), converged = logical(), n = integer())
  out <- new_effect_table(study = rep(info$id, nrow(tab)), feature = tab$feature, term = tab$term,
                          estimate = tab$estimate, se = tab$se, rank = tab$rank, part = tab$part,
                          df = tab$df, p_value = tab$p_value, n = tab$n, engine = rep(engine, nrow(tab)),
                          effect_scale = tab$effect_scale, converged = tab$converged)
  attr(out, "fit_info") <- fit_info
  attr(out, "low_support") <- low_support
  attr(out, "fit_settings") <- list(formula = deparse(formula), engine = engine, rank = rank, longitudinal = longitudinal,
                                    min_prevalence = min_prevalence, min_mean = min_mean, exclude_unassigned = exclude_unassigned,
                                    n_cycles = if (engine == "bezi_gamlss") n_cycles else NULL,
                                    prevalence = if (!is.null(prevalence)) deparse(prevalence),
                                    zero_inflation = if (engine %in% c("zibb_glmmtmb", "bezi_glmmtmb")) zero_inflation,
                                    pseudocount = pseudocount, min_nonzero = need_nz)
  if (!is.null(fit_info) && any(!is.na(fit_info$error))) {
    warning(sum(!is.na(fit_info$error)), " of ", nrow(fit_info), " fits failed (first error: ", fit_info$error[!is.na(fit_info$error)][1],
            "); see attr(, \"fit_info\").", call. = FALSE)
  }
  if (!is.null(fit_info) && any(!fit_info$converged & is.na(fit_info$error))) {
    nc <- sum(!fit_info$converged & is.na(fit_info$error))
    prev_fit <- stats::setNames(vapply(idx, function(i) mean(ra[i, ] > 0, na.rm = TRUE), numeric(1)), rownames(ra)[idx])
    bad <- !fit_info$converged & is.na(fit_info$error)
    warning(nc, " of ", nrow(fit_info), " fits did not converge",
            if (engine == "bezi_gamlss") paste0(" within ", n_cycles, " iterations; try a larger 'n_cycles'"), "; see attr(, \"fit_info\").",
            .convergence_hint(prev_fit[fit_info$feature[bad]], prev_fit[fit_info$feature[fit_info$converged]], min_prevalence), call. = FALSE)
  }
  if (!is.null(fit_info) && any(fit_info$singular)) {
    .inform(sum(fit_info$singular), " of ", nrow(fit_info), " mixed-model fits were singular (the subject variance is estimated as zero), ",
            "so those estimates equal a model without the subject effect.")
  }
  out
}

# Centred log-ratio matrix (features x samples) over the features that make up the composition
.clr_matrix <- function(study, composition, pseudocount) {
  counts <- mm_counts(study); ra <- mm_relabund(study)
  X <- if (!is.null(counts)) counts[composition, , drop = FALSE] + (if (is.null(pseudocount)) 0.5 else pseudocount) else {
    r <- ra[composition, , drop = FALSE]
    r + (if (is.null(pseudocount)) min(r[r > 0], na.rm = TRUE) / 2 else pseudocount)
  }
  L <- log(X)
  out <- sweep(L, 2, colMeans(L), "-")
  out[, mm_depth(study) %in% 0 | colSums(is.na(ra)) > 0] <- NA_real_
  out
}

.coef_rows <- function(ct, part, scale, converged, p_col = "Pr(>|t|)") {
  ct <- ct[rownames(ct) != "(Intercept)", , drop = FALSE]
  if (!nrow(ct)) return(NULL)
  data.frame(term = rownames(ct), part = part, estimate = unname(ct[, "Estimate"]), se = unname(ct[, "Std. Error"]),
             p_value = unname(ct[, p_col]), df = NA_real_, effect_scale = scale, converged = converged,
             stringsAsFactors = FALSE)
}

# One GAMLSS-BEZI fit; mean-part coefficients without the intercept
.fit_bezi_gamlss <- function(f, d, n_cycles = 100) {
  warns <- character(); msgs <- character()
  catch <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) {
    msgs <<- c(msgs, trimws(conditionMessage(m))); invokeRestart("muffleMessage")
  })
  fit <- tryCatch(catch(gamlss::gamlss(f, family = "BEZI", data = d, trace = FALSE,
                                       control = gamlss::gamlss.control(n.cyc = n_cycles, trace = FALSE))),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(rows = NULL, converged = FALSE, warnings = c(warns, msgs[!grepl("singular", msgs)]), singular = any(grepl("singular", msgs)), error = conditionMessage(fit)))
  }
  sm <- NULL
  utils::capture.output(sm <- tryCatch(catch(summary(fit, save = TRUE)), error = function(e) e))
  if (inherits(sm, "error")) {
    return(list(rows = NULL, converged = isTRUE(fit$converged), warnings = c(warns, msgs[!grepl("singular", msgs)]), singular = any(grepl("singular", msgs)), error = conditionMessage(sm)))
  }
  # The table stacks the mean, sigma and nu parts, each starting with an intercept. Aliased mean coefficients are not
  # in it, so count the rows up to the second intercept and not the coefficients of the fit.
  ints <- which(rownames(sm$coef.table) == "(Intercept)")
  k <- if (length(ints) >= 2) ints[2] - 1L else nrow(sm$coef.table)
  conv <- isTRUE(fit$converged)
  list(rows = .coef_rows(sm$coef.table[seq_len(k), , drop = FALSE], "abundance", .mm_engine_scale[["bezi_gamlss"]], conv),
       converged = conv, warnings = c(warns, msgs[!grepl("singular", msgs)]), singular = any(grepl("singular", msgs)), error = NA_character_)
}

# One linear (mixed) model on a transformed scale
.fit_lm <- function(f, d, longitudinal, scale) {
  warns <- character(); msgs <- character()
  catch <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) {
    msgs <<- c(msgs, trimws(conditionMessage(m))); invokeRestart("muffleMessage")
  })
  fit <- tryCatch(
    if (longitudinal) catch(lme4::lmer(f, data = d)) else catch(stats::lm(f, data = d)),
    error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(rows = NULL, converged = FALSE, warnings = c(warns, msgs[!grepl("singular", msgs)]), singular = any(grepl("singular", msgs)), error = conditionMessage(fit)))
  }
  ct <- catch(summary(fit))$coefficients
  if (longitudinal) {
    # lme4 gives no p-values: normal approximation
    ct <- cbind(ct[, c("Estimate", "Std. Error", "t value"), drop = FALSE],
                `Pr(>|t|)` = 2 * stats::pnorm(-abs(ct[, "t value"])))
    conv <- isTRUE(fit@optinfo$conv$opt == 0) && !any(grepl("converge", warns))
    df <- NA_real_
  } else {
    conv <- TRUE
    df <- as.numeric(stats::df.residual(fit))
  }
  rows <- .coef_rows(ct, "abundance", scale, conv)
  if (!is.null(rows)) rows$df <- df
  list(rows = rows, converged = conv, warnings = c(warns, msgs[!grepl("singular", msgs)]), singular = any(grepl("singular", msgs)), error = NA_character_)
}

# One zero-inflated beta-binomial fit (glmmTMB)
.fit_zibb <- function(terms_f, d, longitudinal, prevalence, zero_inflation) {
  warns <- character(); msgs <- character()
  catch <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) {
    msgs <<- c(msgs, trimws(conditionMessage(m))); invokeRestart("muffleMessage")
  })
  use_zi <- isTRUE(zero_inflation) && any(d$.y == 0)
  f <- stats::reformulate(c(terms_f, if (longitudinal) "(1 | .subject)"), response = "cbind(.y, .fail)")
  zif <- if (!use_zi) ~0 else if (is.null(prevalence)) ~1 else prevalence
  fit <- tryCatch(catch(glmmTMB::glmmTMB(f, ziformula = zif, family = glmmTMB::betabinomial(), data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(rows = NULL, converged = FALSE, warnings = c(warns, msgs[!grepl("singular", msgs)]), singular = any(grepl("singular", msgs)), error = conditionMessage(fit)))
  }
  conv <- isTRUE(fit$fit$convergence == 0) && isTRUE(fit$sdr$pdHess)
  sm <- tryCatch(catch(summary(fit))$coefficients, error = function(e) e)
  if (inherits(sm, "error")) {
    return(list(rows = NULL, converged = conv, warnings = c(warns, msgs[!grepl("singular", msgs)]), singular = any(grepl("singular", msgs)), error = conditionMessage(sm)))
  }
  rows <- .coef_rows(sm$cond, "abundance", .mm_engine_scale[["zibb_glmmtmb"]], conv, "Pr(>|z|)")
  if (use_zi && !is.null(prevalence) && !is.null(sm$zi)) {
    zr <- .coef_rows(sm$zi, "prevalence", "logit_zero_prob", conv, "Pr(>|z|)")
    if (!is.null(zr)) {
      zr$converged <- conv & is.finite(zr$se) & zr$se <= .zi_se_limit
      rows <- rbind(rows, zr)
    }
  }
  list(rows = rows, converged = conv, warnings = c(warns, msgs[!grepl("singular", msgs)]), singular = any(grepl("singular", msgs)), error = NA_character_)
}

# One zero-inflated beta fit on relative abundance (glmmTMB)
.fit_bezi_glmmtmb <- function(terms_f, d, longitudinal, prevalence, zero_inflation) {
  warns <- character(); msgs <- character()
  catch <- function(expr) withCallingHandlers(expr, warning = function(w) {
    warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) {
    msgs <<- c(msgs, trimws(conditionMessage(m))); invokeRestart("muffleMessage")
  })
  has_zero <- any(d$.y == 0)
  use_zi <- isTRUE(zero_inflation) && has_zero
  if (has_zero && !use_zi) d$.y[d$.y == 0] <- min(d$.y[d$.y > 0]) / 2       # a plain beta model needs values above zero
  f <- stats::reformulate(c(terms_f, if (longitudinal) "(1 | .subject)"), response = ".y")
  zif <- if (!use_zi) ~0 else if (is.null(prevalence)) ~1 else prevalence
  fit <- tryCatch(catch(glmmTMB::glmmTMB(f, ziformula = zif, family = glmmTMB::beta_family(), data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    return(list(rows = NULL, converged = FALSE, warnings = c(warns, msgs[!grepl("singular", msgs)]), singular = any(grepl("singular", msgs)), error = conditionMessage(fit)))
  }
  conv <- isTRUE(fit$fit$convergence == 0) && isTRUE(fit$sdr$pdHess)
  sing <- longitudinal && isTRUE(tryCatch(unname(glmmTMB::VarCorr(fit)$cond[[1]][1, 1]) < 1e-8, error = function(e) FALSE))
  sm <- tryCatch(catch(summary(fit))$coefficients, error = function(e) e)
  if (inherits(sm, "error")) {
    return(list(rows = NULL, converged = conv, warnings = c(warns, msgs[!grepl("singular", msgs)]), singular = sing, error = conditionMessage(sm)))
  }
  rows <- .coef_rows(sm$cond, "abundance", .mm_engine_scale[["bezi_glmmtmb"]], conv, "Pr(>|z|)")
  if (use_zi && !is.null(prevalence) && !is.null(sm$zi)) {
    zr <- .coef_rows(sm$zi, "prevalence", "logit_zero_prob", conv, "Pr(>|z|)")
    if (!is.null(zr)) {
      zr$converged <- conv & is.finite(zr$se) & zr$se <= .zi_se_limit
      rows <- rbind(rows, zr)
    }
  }
  list(rows = rows, converged = conv, warnings = c(warns, msgs[!grepl("singular", msgs)]), singular = sing, error = NA_character_)
}

# All features at once: they share the design matrix, so a single QR factorisation gives every coefficient, standard
# error and p-value exactly as lm() would. Returns NULL (use the per-feature fits) if the design is not full rank, and
# marks features with a (nearly) perfect fit, for which lm() warns, so that they go through lm() as before.
.lm_batch <- function(rhs, d0, Yb) {
  if (anyNA(Yb)) return(NULL)
  X <- tryCatch(stats::model.matrix(rhs, d0), error = function(e) NULL)       # e.g. a variable with one level: let the per-feature fits report it
  if (is.null(X)) return(NULL)
  qrX <- qr(X)
  rdf <- nrow(X) - ncol(X)
  if (qrX$rank < ncol(X) || rdf < 1 || !identical(as.integer(qrX$pivot), seq_len(ncol(X)))) return(NULL)
  B <- qr.coef(qrX, t(Yb)); R <- t(Yb) - X %*% B
  s2 <- colSums(R^2) / rdf
  se <- sqrt(outer(diag(chol2inv(qr.R(qrX))), s2)); dimnames(se) <- dimnames(B); tv <- B / se
  list(B = B, se = se, tv = tv, pv = 2 * stats::pt(-abs(tv), rdf), rdf = as.numeric(rdf),
       ok = is.finite(s2) & s2 > 1e-20 * pmax(1, rowMeans(Yb^2)))
}

# When fits fail, they usually fail for the rare features; say so, with the user's own numbers
.convergence_hint <- function(prev_failed, prev_ok, min_prevalence) {
  prev_failed <- prev_failed[is.finite(prev_failed)]; prev_ok <- prev_ok[is.finite(prev_ok)]
  if (length(prev_failed) < 5 || !length(prev_ok) || stats::median(prev_failed) >= stats::median(prev_ok)) return("")
  suggest <- min(0.5, max(2 * min_prevalence, ceiling(stats::quantile(prev_failed, 0.75, names = FALSE) / 0.05) * 0.05))
  paste0(" The features that did not converge are rarer (present in a median of ", .pct(stats::median(prev_failed)), " of the samples, against ",
         .pct(stats::median(prev_ok)), " for the others); try a higher 'min_prevalence', for example ", format(suggest), ".")
}
