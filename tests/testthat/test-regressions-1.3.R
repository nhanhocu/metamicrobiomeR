# Regression tests for the fixes in metamicrobiomeR 1.3

# Small example from the Bangladesh data (5 fusobacteria taxa)
example_data <- function() {
  utils::data("taxtab6", package = "metamicrobiomeR", envir = environment())
  # as.data.frame() keeps this helper usable on the old grouped_df format too
  taxtab6 <- as.data.frame(get("taxtab6", envir = environment()))
  tl <- colnames(taxtab6)[grep("k__bacteria.p__fusobacteria", colnames(taxtab6))]
  taxtab6[, c("personid", "x.sampleid", "bf", "age.sample", tl)]
}

fit_lmm <- function(...) {
  suppressMessages(suppressWarnings(
    taxa.compare(taxtab = example_data(), propmed.rel = "lm", comvar = "bf",
                 adjustvar = "age.sample", longitudinal = "yes",
                 p.adjust.method = "fdr", ...)
  ))
}

test_that("bundled datasets are plain data frames that subset under dplyr", {
  for (nm in c("taxtab6", "asum4", "covar.rm", "tabsex4")) {
    utils::data(list = nm, package = "metamicrobiomeR", envir = environment())
    x <- get(nm, envir = environment())
    expect_false(inherits(x, "grouped_df"), info = nm)
    expect_s3_class(x, "data.frame")
    expect_no_error(x[, 1:2])
  }
})

test_that("mixed-model p-values are two-sided normal-approximation p-values", {
  skip_on_cran()
  res <- fit_lmm(transform = "none")
  est <- res[["Estimate.age.sample"]]
  se <- res[["Std. Error.age.sample"]]
  expect_equal(res[["Pr(>|t|).age.sample"]], 2 * stats::pnorm(-abs(est / se)))
})

test_that("adjusted p-value columns keep the full covariate name", {
  skip_on_cran()
  res <- fit_lmm(transform = "none")
  expect_true("pval.adjust.age.sample" %in% colnames(res))
  expect_false("pval.adjust.sample" %in% colnames(res))
})

test_that("GMPR transform runs and gives finite size-factor-scaled output", {
  skip_on_cran()
  res <- fit_lmm(transform = "gmpr")
  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0)
  expect_true(any(is.finite(res[["Estimate.age.sample"]])))
})
