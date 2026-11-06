# Tests for R/90-simulate-data.R, R/81-report-manifest.R
test_that("simulate_study returns a study with the requested size and known truth", {
  s <- simulate_study(n_samples = 50, n_features = 20, n_diff = 4, effect = 1.2, seed = 1)
  expect_s3_class(s$study, "mm_study"); expect_equal(dim(mm_counts(s$study)), c(20, 50))
  expect_equal(sum(s$truth$diff), 4); expect_equal(sum(s$truth$beta > 0), 2)            # balanced: alternating directions
  expect_equal(unique(abs(s$truth$beta[s$truth$diff])), 1.2)
  u <- simulate_study(n_features = 20, n_diff = 4, balanced = FALSE, seed = 1); expect_true(all(u$truth$beta[u$truth$diff] > 0))
  expect_equal(simulate_study(seed = 5)$study$assays$counts, simulate_study(seed = 5)$study$assays$counts)
  expect_error(simulate_study(n_features = 3, n_diff = 5), "cannot exceed")
  d <- mm_depth(simulate_study(n_samples = 200, depth = 10000, depth_sd = 0.6, seed = 2)$study)
  expect_gt(max(d) / min(d), 3)                                                           # depth varies
  z <- simulate_study(zero_inflation = 0.3, seed = 3)$study; z0 <- simulate_study(zero_inflation = 0, seed = 3)$study
  expect_gt(mean(mm_counts(z) == 0), mean(mm_counts(z0) == 0) + 0.15)
})

test_that("the simulated group effect is recoverable and compositional effects appear when unbalanced", {
  s <- simulate_study(n_samples = 200, n_features = 12, n_diff = 4, effect = 1.5, depth_sd = 0.3, seed = 4)
  e <- fit_study(s$study, ~ group, engine = "lm_clr", longitudinal = FALSE, min_prevalence = 0)
  est <- e$estimate[e$term == "groupb"]; names(est) <- e$feature[e$term == "groupb"]
  di <- s$truth$feature[s$truth$diff]
  expect_equal(unname(sign(est[di])), sign(s$truth$beta[s$truth$diff]))
  expect_gt(mean(abs(est[di])), 2 * mean(abs(est[setdiff(names(est), di)])))
  # unbalanced: the non-differential features shift in the opposite direction
  u <- simulate_study(n_samples = 200, n_features = 12, n_diff = 6, effect = 2, balanced = FALSE, depth_sd = 0.3, seed = 4)
  r <- mm_relabund(u$study); g <- mm_samples(u$study)$group
  nd <- u$truth$feature[!u$truth$diff]
  expect_lt(mean(rowMeans(r[nd, g == "b"]) / rowMeans(r[nd, g == "a"])), 0.9)
})

test_that("repeated measures and several studies are simulated consistently", {
  s <- simulate_study(n_samples = 40, n_subjects = 20, between = TRUE, seed = 1)$study
  expect_equal(mm_study_info(s)$subject, "person")
  expect_true(all(tapply(mm_samples(s)$group, mm_samples(s)$person, function(x) length(unique(x))) == 1))
  m <- simulate_studies(n_studies = 3, tau = 0.5, small_study = TRUE, n_samples = 50, n_features = 15, seed = 2)
  expect_length(m$studies, 3); expect_equal(sapply(m$studies, function(x) ncol(mm_counts(x))), c(50, 50, 10))
  expect_equal(vapply(m$studies, function(x) mm_study_info(x)$id, ""), c("sim1", "sim2", "sim3"))
  expect_equal(sum(m$truth$diff), 6)
  eff <- lapply(m$studies, function(x) fit_study(x, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0, min_nonzero = 0))
  expect_s3_class(pool_effects(bind_effects(eff), term = "groupb"), "mm_pooled")
})

test_that("manifests describe studies, effects and pooled results", {
  s <- simulate_study(n_samples = 30, n_features = 10, seed = 1)$study
  m <- manifest(s)
  expect_s3_class(m, "mm_manifest"); expect_equal(m$features, 10); expect_equal(m$samples, 30); expect_equal(m$study$id, "sim")
  expect_true(all(c("min", "median", "max") %in% names(m$depth))); expect_equal(m$zero_depth_samples, 0)
  expect_true(!is.na(m$environment$R))
  e <- fit_study(s, ~ group, engine = "lm_asin", longitudinal = FALSE)
  me <- manifest(e)
  expect_equal(me$engines, "lm_asin"); expect_equal(me$fit_settings$formula, "~group"); expect_equal(me$fit_settings$engine, "lm_asin")
  e2 <- e; e2$study <- "other"; class(e2) <- class(e)
  p <- pool_effects(bind_effects(e, e2), method = "DL", test = "z")
  mp <- manifest(p); expect_equal(mp$settings$method, "DL"); expect_equal(mp$studies_per_feature[["max"]], 2)
  f <- tempfile(fileext = ".txt"); write_manifest(s, f)
  txt <- readLines(f); expect_true(any(grepl("^features: 10", txt))); expect_true(any(grepl("environment:", txt)))
  expect_output(print(m), "samples: 30")
  skip_if_not_installed("jsonlite")
  j <- tempfile(fileext = ".json"); write_manifest(e, j, "json")
  expect_equal(jsonlite::fromJSON(j)$type, "effects")
  expect_error(manifest(1), "mm_study")
})
