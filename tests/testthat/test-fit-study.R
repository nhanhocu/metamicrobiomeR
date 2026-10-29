make_bezi_study <- function(n = 90, seed = 11, subject = NULL) {
  set.seed(seed)
  grp <- factor(rep(c("a", "b"), length.out = n))
  age <- runif(n, 0, 6)
  mu1 <- stats::plogis(-1.5 + 1.2 * (grp == "b"))        # feature with a group effect
  mu2 <- rep(stats::plogis(-2), n)                       # feature without
  y1 <- gamlss.dist::rBEZI(n, mu = mu1, sigma = 1, nu = 0.05)
  y2 <- gamlss.dist::rBEZI(n, mu = mu2, sigma = 1, nu = 0.05)
  rare <- c(rep(0, n - 2), 0.3, 0.2)                     # present in 2 samples only
  rel <- rbind(F1 = y1, F2 = y2, Rare = rare) * 0.3
  colnames(rel) <- paste0("S", seq_len(n))
  smp <- data.frame(group = grp, age = age, row.names = colnames(rel))
  if (!is.null(subject)) smp$person <- rep(seq_len(n / 2), each = 2)
  as_study(relabund = rel, samples = smp,
           study = list(id = "sim", feature_type = "taxon", subject = subject))
}

test_that("fit_study returns an effect table with the documented content", {
  skip_if_not_installed("gamlss")
  s <- make_bezi_study()
  e <- suppressWarnings(fit_study(s, ~ group + age))
  expect_s3_class(e, "mm_effects")
  expect_setequal(unique(e$term), c("groupb", "age"))
  expect_true(all(e$part == "abundance" & e$engine == "bezi_gamlss" & e$effect_scale == "logit_mean"))
  expect_true(all(e$study == "sim"))
  g1 <- e[e$feature == "F1" & e$term == "groupb", ]
  expect_gt(g1$estimate, 0.3)
  expect_lt(g1$p_value, 0.05)
  expect_gt(g1$se, 0)
  expect_equal(g1$n, 90L)
  info <- attr(e, "fit_info")
  expect_setequal(info$feature, c("F1", "F2"))
})

test_that("features failing the prevalence filter are not fitted", {
  skip_if_not_installed("gamlss")
  e <- suppressWarnings(fit_study(make_bezi_study(), ~ group))
  expect_false("Rare" %in% e$feature)
  e2 <- suppressWarnings(fit_study(make_bezi_study(), ~ group, min_prevalence = 0.01))
  expect_true("Rare" %in% attr(e2, "fit_info")$feature)
})

test_that("a longitudinal model adds a subject random intercept", {
  skip_if_not_installed("gamlss")
  s <- make_bezi_study(subject = "person")
  e <- suppressWarnings(fit_study(s, ~ group))
  expect_type(e$converged, "logical")
  # the feature with a real group effect converges; convergence of the null
  # feature is not guaranteed with this engine (see dev/VALIDATION.md)
  expect_true(e$converged[e$feature == "F1"])
  expect_setequal(unique(e$term), "groupb")
  expect_gt(e$estimate[e$feature == "F1"], 0.3)
  expect_error(fit_study(make_bezi_study(), ~ group, longitudinal = TRUE), "subject")
})

test_that("fit_study validates its input", {
  s <- make_bezi_study()
  expect_error(fit_study(s, group ~ age), "one-sided")
  expect_error(fit_study(s, ~ nope), "not in the sample data")
  expect_error(fit_study(s, ~ group, rank = "genus"), "rank")
  expect_error(fit_study(list(), ~ group), "mm_study")
})

legacy_ref <- function() readRDS(testthat::test_path("fixtures", "legacy-reference.rds"))
legacy_study <- function(ref) {
  tab <- ref$tab; phyla <- grep("^k__", colnames(tab), value = TRUE)
  rel <- t(as.matrix(tab[, phyla])); colnames(rel) <- tab$x.sampleid
  smp <- tab[, c("bf", "age.sample", "personid")]; rownames(smp) <- colnames(rel)
  as_study(relabund = rel, samples = smp, study = list(id = "bangladesh", feature_type = "otu", subject = "personid"))
}
expect_matches_legacy <- function(new, old) {
  for (term in c("bfNon_exclusiveBF", "bfNo_BF", "age.sample")) {
    n <- new[new$term == term, ]; n <- n[match(old$id, n$feature), ]
    expect_equal(n$estimate, old[[paste0("Estimate.", term)]], tolerance = 1e-8)
    expect_equal(n$se, old[[paste0("Std. Error.", term)]], tolerance = 1e-8)
    expect_equal(n$p_value, old[[paste0("Pr(>|t|).", term)]], tolerance = 1e-8)
  }
}

test_that("fit_study reproduces the 1.x GAMLSS-BEZI results (stored reference, Bangladesh data)", {
  skip_if_not_installed("gamlss")
  ref <- legacy_ref()
  new <- suppressWarnings(fit_study(legacy_study(ref), ~ bf + age.sample, n_cycles = 20))
  expect_equal(length(unique(new$feature)), 4)
  expect_matches_legacy(new, ref$bezi)
})

test_that("a warning reports fits that did not converge", {
  skip_if_not_installed("gamlss")
  s <- make_bezi_study()
  expect_warning(fit_study(s, ~ group, n_cycles = 2), "did not converge")
})

test_that("features unassigned at their own rank are excluded unless asked for", {
  skip_if_not_installed("gamlss")
  s <- make_bezi_study()
  tx <- data.frame(genus = c("Gen1", NA, "Gen3"), rank = "genus", row.names = rownames(mm_relabund(s)))
  s2 <- as_study(relabund = mm_relabund(s), taxonomy = tx, samples = mm_samples(s),
                 study = list(id = "sim", feature_type = "taxon"))
  e <- suppressWarnings(fit_study(s2, ~ group))
  expect_false("F2" %in% e$feature)
  e2 <- suppressWarnings(fit_study(s2, ~ group, exclude_unassigned = FALSE))
  expect_true("F2" %in% e2$feature)
})

test_that("the lm_asin engine returns effects on the arcsine square root scale", {
  s <- make_bezi_study()
  e <- fit_study(s, ~ group + age, engine = "lm_asin")
  expect_true(all(e$engine == "lm_asin" & e$effect_scale == "asin_sqrt_mean"))
  expect_setequal(unique(e$term), c("groupb", "age"))
  g1 <- e[e$feature == "F1" & e$term == "groupb", ]
  expect_gt(g1$estimate, 0.05)
  expect_lt(g1$p_value, 0.05)
  expect_gt(g1$df, 0)                       # cross-sectional: lm with t-based p-values
  s2 <- make_bezi_study(subject = "person")
  e2 <- fit_study(s2, ~ group, engine = "lm_asin")
  expect_true(all(is.na(e2$df)))            # mixed model: normal approximation, no df
  expect_equal(e2$p_value, 2 * stats::pnorm(-abs(e2$estimate / e2$se)))
  other <- suppressWarnings(fit_study(s, ~ group + age)); other$study <- "sim2"   # a different study, bezi engine
  expect_error(bind_effects(e, other), "effect scales")
})

test_that("lm_asin reproduces the 1.x linear mixed and linear models (stored reference)", {
  ref <- legacy_ref(); s <- legacy_study(ref)
  expect_matches_legacy(fit_study(s, ~ bf + age.sample, engine = "lm_asin", longitudinal = TRUE), ref$lmasin_mixed)
  expect_matches_legacy(fit_study(s, ~ bf + age.sample, engine = "lm_asin", longitudinal = FALSE), ref$lmasin_linear)
})
