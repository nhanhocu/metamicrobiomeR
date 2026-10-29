mkeff <- function(study, est, se, feature = "g1", term = "x", part = "abundance") {
  new_effect_table(study = study, feature = feature, term = term, estimate = est, se = se, rank = "genus", part = part,
                   n = 30L, engine = "t", effect_scale = "logit_mean")
}
studies <- paste0("s", 1:6)
sd_ <- data.frame(study = studies, region = c("V4", "V4", "V4", "V1V2", "V1V2", "V1V2"),
                  year = c(2012, 2014, 2016, 2013, 2015, 2017), stringsAsFactors = FALSE)
est <- c(0.20, 0.25, 0.15, 0.80, 0.70, 0.85); se <- c(0.10, 0.12, 0.10, 0.15, 0.10, 0.12)
eff <- bind_effects(mkeff(studies, est, se))

test_that("meta_regress matches a direct metafor fit and reports the moderator effect", {
  m <- meta_regress(eff, sd_, ~ region, method = "REML", test = "z")
  expect_s3_class(m, "mm_metareg")
  d <- cbind(data.frame(yi = est, sei = se), sd_)
  fit <- metafor::rma.uni(yi, sei = sei, mods = ~ region, data = d, method = "REML", test = "z")
  expect_equal(m$moderator, "regionV4")
  expect_equal(m$estimate, as.numeric(fit$beta)[2], tolerance = 1e-8)
  expect_equal(m$se, fit$se[2], tolerance = 1e-8)
  expect_equal(m$p_value, fit$pval[2], tolerance = 1e-8)
  expect_equal(m$QM, fit$QM, tolerance = 1e-8)
  expect_lt(m$estimate, -0.4)                                   # V4 studies have the smaller effect
  expect_equal(m$k, 6)
})

test_that("a numeric moderator gives one slope, and several moderators give several rows", {
  expect_equal(meta_regress(eff, sd_, ~ year, test = "z")$moderator, "year")
  m2 <- meta_regress(eff, sd_, ~ region + year, test = "z")
  expect_equal(nrow(m2), 2)
  expect_equal(length(unique(m2$QM)), 1)                        # one omnibus test per feature
})

test_that("features with too few studies or a moderator without variation are skipped", {
  few <- bind_effects(mkeff(studies, est, se, "g1"), mkeff(studies[1:2], c(0.1, 0.2), c(0.1, 0.1), "g2"))
  m <- meta_regress(few, sd_, ~ region, test = "z")
  expect_equal(unique(m$feature), "g1")                         # g2 has 2 studies, both V4
  expect_error(meta_regress(bind_effects(mkeff(studies[1:2], c(0.1, 0.2), c(0.1, 0.1))), sd_, ~ region), "enough studies")
})

test_that("input is validated", {
  expect_error(meta_regress(data.frame(a = 1), sd_, ~ region), "mm_effects")
  expect_error(meta_regress(eff, sd_[, "region", drop = FALSE], ~ region), "column 'study'")
  expect_error(meta_regress(eff, sd_, ~ nope), "not in 'study_data'")
  expect_error(meta_regress(eff, sd_[1:3, ], ~ region), "No study_data")
  expect_error(meta_regress(eff, sd_, "region"), "one-sided")
})

test_that("the omnibus p-value is adjusted across features and study_table summarises studies", {
  e2 <- bind_effects(mkeff(studies, est, se, "g1"), mkeff(studies, est[c(2, 1, 3, 5, 4, 6)], se, "g2"), mkeff(studies, rev(est), se, "g3"))
  m <- meta_regress(e2, sd_, ~ region, test = "z")
  one <- m[!duplicated(m$feature), ]
  expect_equal(sort(one$QM_p_adjust), sort(stats::p.adjust(one$QM_p, "fdr")), tolerance = 1e-10)
  s1 <- as_study(counts = matrix(1:6, 2, 3, dimnames = list(c("a", "b"), c("x", "y", "z"))), samples = data.frame(g = 1:3, row.names = c("x", "y", "z")),
                 study = list(id = "A", feature_type = "asv", region = "V4", platform = "MiSeq"))
  s2 <- as_study(counts = matrix(1:4, 2, 2, dimnames = list(c("a", "b"), c("x", "y"))), samples = data.frame(g = 1:2, row.names = c("x", "y")),
                 study = list(id = "B", feature_type = "asv", region = "V1V2", measured = FALSE))
  st <- study_table(s1, s2)
  expect_equal(st$study, c("A", "B")); expect_equal(st$region, c("V4", "V1V2")); expect_equal(st$n_samples, c(3, 2))
  expect_identical(st$measured[2], FALSE)
  expect_error(study_table(s1, s1), "unique")
})
