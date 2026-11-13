# Tests for the example data (R/90-simulate-data.R, data/example_studies.rda) and the claims made in the vignettes
planted <- c("Escherichia", "Enterococcus", "Faecalibacterium", "Roseburia", "Butyricicoccus")

test_that("example_studies is a valid, realistic cross-platform data set", {
  expect_type(example_studies, "list"); expect_named(example_studies, c("CohortA", "CohortB", "CohortC", "CohortD"))
  for (s in example_studies) { expect_s3_class(s, "mm_study"); expect_silent(check_study(s)) }
  st <- study_table(example_studies)
  expect_equal(st$n_samples, rep(50, 4)); expect_equal(length(unique(st$platform)), 3); expect_equal(length(unique(st$database)), 2)
  expect_equal(unname(unlist(lapply(example_studies, function(s) mm_study_info(s)$feature_type))), c("asv", "asv", "taxon", "asv"))
  for (s in example_studies) expect_true(all(c("group", "age", "sex") %in% colnames(mm_samples(s))))
  expect_equal(levels(mm_samples(example_studies$CohortA)$group), c("control", "case"))
  pc <- check_poolable(example_studies, rank = "genus")
  expect_true(all(c("taxonomy databases differ", "platforms differ") %in% pc$issues$issue)); expect_true(pc$ok)
  cov <- feature_coverage(example_studies, "genus")
  expect_equal(cov$k[cov$feature == "Fusobacterium"], 1)                                          # in one study only: cannot be pooled
  expect_gte(cov$k[cov$feature == "Eubacterium"], 3)                                              # "[Eubacterium]" is matched with "Eubacterium"
  expect_gte(cov$k[cov$feature == "Blautia"], 3)                                                  # "Blautia_A" is matched with "Blautia"
})

test_that("the planted differences are what the documentation says they are", {
  res <- suppressMessages(meta_analyse(example_studies, ~ group + age, term = "groupcase", rank = "genus", engine = "lm_asin",
                                       fit_args = list(longitudinal = FALSE), quality = "off", check = FALSE))
  p <- as.data.frame(res); p <- p[order(p$p_value), ]
  expect_equal(nrow(p), 57); expect_false("Fusobacterium" %in% p$feature)
  expect_setequal(p$feature[1:4], c("Faecalibacterium", "Butyricicoccus", "Enterococcus", "Roseburia"))   # the vignette says the four strongest are planted
  est <- stats::setNames(p$estimate, p$feature)
  expect_true(all(est[c("Escherichia", "Enterococcus")] > 0)); expect_true(all(est[c("Faecalibacterium", "Roseburia", "Butyricicoccus")] < 0))
  expect_lt(stats::median(p$I2), 25)
  expect_equal(sum(p$p_adjust <= 0.05), 0)                                                         # the default test is conservative here, as the vignette says
  s <- pool_sensitivity(res$effects, term = "groupcase", significance = "p_adjust")
  expect_equal(unname(s$n_significant), c(47, 47, 13, 0))                                          # the numbers quoted in the vignette
})

test_that("the cross-platform vignette's claims about two-stage and combined analyses hold", {
  res <- suppressMessages(meta_analyse(example_studies, ~ group + age, term = "groupcase", rank = "genus", engine = "lm_asin",
                                       fit_args = list(longitudinal = FALSE), quality = "off", check = FALSE))
  two <- as.data.frame(res)
  merged <- suppressMessages(merge_studies(example_studies, rank = "genus"))
  one <- fit_study(merged, ~ group + age + study, engine = "lm_asin", longitudinal = FALSE); one <- one[one$term == "groupcase", ]
  expect_equal(c(nrow(two), nrow(one)), c(57, 48))
  both <- merge(two[, c("feature", "estimate")], one[, c("feature", "estimate")], by = "feature", suffixes = c("_m", "_c"))
  expect_equal(nrow(both), 48); expect_gt(stats::cor(both$estimate_m, both$estimate_c), 0.99); expect_equal(mean(sign(both$estimate_m) == sign(both$estimate_c)), 1)
  expect_setequal(setdiff(two$feature, one$feature), c("Sutterella", "Lachnospira", "Bacteroides", "Dorea", "Coprococcus", "Bilophila", "Desulfovibrio", "Akkermansia", "Haemophilus"))
  planted <- c("Escherichia", "Enterococcus", "Faecalibacterium", "Roseburia", "Butyricicoccus")
  expect_equal(sum(head(one$feature[order(one$p_value)], 5) %in% planted), 4)
  expect_equal(sum(p.adjust(one$p_value, "fdr") <= 0.05), 43); expect_equal(sum(two$p_adjust <= 0.05), 0)
})
