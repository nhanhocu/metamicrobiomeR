# Tests for R/50-pool-effects.R, R/52-pool-meta-regression.R, R/53-pool-sensitivity.R, R/54-pool-hierarchical.R

# ================ pooling ================
eff <- function(study, est, se, feature = "g1", term = "x", scale = "logit_mean", rank = "genus", part = "abundance") {
  new_effect_table(study = study, feature = feature, term = term, estimate = est, se = se,
                   rank = rank, part = part, n = 30L, engine = "test", effect_scale = scale)
}

# DerSimonian-Laird computed by hand
dl <- function(y, se) {
  w <- 1 / se^2; k <- length(y)
  ybar <- sum(w * y) / sum(w); Q <- sum(w * (y - ybar)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (se^2 + tau2)
  est <- sum(ws * y) / sum(ws); s <- 1 / sqrt(sum(ws))
  list(est = est, se = s, tau2 = tau2, Q = Q, p = 2 * stats::pnorm(-abs(est / s)),
       lo = est - qnorm(0.975) * s, hi = est + qnorm(0.975) * s)
}
y <- c(0.8, 0.1, 0.5, -0.2); se <- c(0.2, 0.3, 0.25, 0.4)
bound <- bind_effects(eff("s1", y[1], se[1]), eff("s2", y[2], se[2]), eff("s3", y[3], se[3]), eff("s4", y[4], se[4]))

test_that("DerSimonian-Laird pooling matches the hand calculation", {
  p <- pool_effects(bound, method = "DL", test = "z")
  h <- dl(y, se)
  expect_s3_class(p, "mm_pooled")
  expect_equal(p$k, 4)
  expect_equal(p$estimate, h$est, tolerance = 1e-8)
  expect_equal(p$se, h$se, tolerance = 1e-8)
  expect_equal(p$tau2, h$tau2, tolerance = 1e-8)
  expect_equal(p$Q, h$Q, tolerance = 1e-8)
  expect_equal(p$p_value, h$p, tolerance = 1e-8)
  expect_equal(c(p$ci_lo, p$ci_hi), c(h$lo, h$hi), tolerance = 1e-8)
})

test_that("the Hartung-Knapp adjustment widens the interval and changes the test", {
  a <- pool_effects(bound, method = "DL", test = "z")
  b <- pool_effects(bound, method = "DL", test = "adhoc")
  expect_equal(a$estimate, b$estimate, tolerance = 1e-8)
  expect_gte(b$ci_hi - b$ci_lo, a$ci_hi - a$ci_lo)
  expect_false(is.na(b$pi_lo))              # prediction interval with 4 studies
  two <- pool_effects(bind_effects(eff("s1", 0.5, 0.2), eff("s2", 0.4, 0.2)), method = "DL")
  expect_true(is.na(two$pi_lo))             # not reported with fewer than 3 studies
})

test_that("plain Hartung-Knapp is anti-conservative with two near-identical studies; the ad hoc default is not", {
  two <- bind_effects(eff("s1", 0.050, 0.5), eff("s2", 0.052, 0.5))     # not significant by any sensible test
  z <- pool_effects(two, method = "REML", test = "z")
  knha <- pool_effects(two, method = "REML", test = "knha")
  adhoc <- pool_effects(two, method = "REML", test = "adhoc")
  expect_lt(knha$p_value, 0.05)             # the failure mode: spuriously small p (z gives p > 0.5)
  expect_gt(z$p_value, 0.5)
  expect_gte(adhoc$p_value, z$p_value)      # the default never beats the normal-theory test
  expect_identical(attr(pool_effects(two), "settings")$test, "adhoc")
})

test_that("features need enough studies, and a missing feature is not pooled", {
  b <- bind_effects(eff("s1", 0.5, 0.2, "g1"), eff("s2", 0.4, 0.2, "g1"), eff("s3", 0.6, 0.2, "g1"),
                    eff("s1", 0.1, 0.2, "g2"),                       # only one study
                    eff("s2", 0.3, 0.2, "g3"), eff("s3", 0.2, 0.2, "g3"))
  p <- pool_effects(b, method = "DL", test = "z")
  expect_setequal(p$feature, c("g1", "g3"))
  expect_equal(p$k[p$feature == "g3"], 2)
  p2 <- pool_effects(b, method = "DL", test = "z", min_prop = 1)       # must be in all 3 studies
  expect_equal(p2$feature, "g1")
})

test_that("p-values are adjusted within the requested family", {
  mk <- function(feature, rank, est) bind_effects(eff("s1", est, 0.1, feature, rank = rank), eff("s2", est, 0.1, feature, rank = rank))
  b <- bind_effects(mk("a", "genus", 0.5), mk("b", "genus", 0.1), mk("c", "family", 0.5), mk("d", "family", 0.1))
  by_rank <- pool_effects(b, method = "DL", test = "z")                                   # default: within rank
  mixed <- pool_effects(b, method = "DL", test = "z", adjust_by = c("part", "term"))
  expect_equal(by_rank$p_adjust[by_rank$feature == "b"], stats::p.adjust(by_rank$p_value[by_rank$rank == "genus"], "fdr")[by_rank$feature[by_rank$rank == "genus"] == "b"])
  expect_equal(sort(mixed$p_adjust), sort(stats::p.adjust(mixed$p_value, "fdr")))
  expect_equal(pool_effects(b, method = "DL", test = "z", p_adjust = "none")$p_adjust,
               pool_effects(b, method = "DL", test = "z", p_adjust = "none")$p_value)
})

test_that("terms and parts are pooled separately and can be selected", {
  b <- bind_effects(eff("s1", 0.5, 0.2, term = "x"), eff("s2", 0.4, 0.2, term = "x"),
                    eff("s1", 0.2, 0.2, term = "z"), eff("s2", 0.1, 0.2, term = "z"))
  expect_setequal(pool_effects(b, method = "DL")$term, c("x", "z"))
  expect_equal(pool_effects(b, term = "z", method = "DL")$term, "z")
  expect_error(pool_effects(b, term = "nope"), "No effects")
})

test_that("invalid input and unusable estimates are handled", {
  expect_error(pool_effects(data.frame(a = 1)), "mm_effects")
  b <- bind_effects(eff("s1", 0.5, 0.2), eff("s2", NA, 0.2), eff("s3", 0.4, 0.2))
  expect_warning(p <- pool_effects(b, method = "DL"), "not used")
  expect_equal(p$k, 2)
  expect_error(pool_effects(bind_effects(eff("s1", 0.5, 0.2)), min_studies = 2), "No feature")
})

test_that("influence diagnostics report the weight of the most influential study", {
  p <- pool_effects(bound, method = "DL", test = "z")
  h <- dl(y, se)
  ws <- 1 / (se^2 + h$tau2); share <- ws / sum(ws)
  expect_equal(p$max_weight, max(share), tolerance = 1e-8)
  expect_equal(p$top_study, "s1")                              # the study with the smallest standard error
  # one very precise study carries almost all the weight
  dom <- bind_effects(eff("a", 0.4, 0.005), eff("b", 0.0, 0.4), eff("c", 0.5, 0.35))
  expect_gt(pool_effects(dom, method = "DL", test = "z")$max_weight, 0.99)
  expect_identical(pool_effects(dom, method = "DL", test = "z")$top_study, "a")
})

test_that("leave_one_out repeats the pooling without each study", {
  mkf <- function(study, feature, est, se) eff(study, est, se, feature = feature)
  b <- bind_effects(mkf("s1", "g", 0.8, 0.2), mkf("s2", "g", 0.1, 0.3), mkf("s3", "g", 0.5, 0.25), mkf("s4", "g", -0.2, 0.4))
  l <- leave_one_out(b, method = "DL", test = "z")
  expect_setequal(l$omitted, c("s1", "s2", "s3", "s4"))
  expect_true(all(l$k == 3))
  for (s in c("s1", "s3")) {
    keep <- setdiff(c("s1", "s2", "s3", "s4"), s)
    h <- dl(c(s1 = 0.8, s2 = 0.1, s3 = 0.5, s4 = -0.2)[keep], c(s1 = 0.2, s2 = 0.3, s3 = 0.25, s4 = 0.4)[keep])
    expect_equal(l$estimate[l$omitted == s], h$est, tolerance = 1e-8)
    expect_equal(l$p_value[l$omitted == s], h$p, tolerance = 1e-8)
  }
  # the result that depends on one study loses significance without it
  dom <- bind_effects(eff("a", 0.43, 0.005), eff("b", 0.01, 0.43), eff("c", 0.55, 0.35))
  l2 <- leave_one_out(dom, method = "DL", test = "z")
  expect_lt(pool_effects(dom, method = "DL", test = "z")$p_value, 1e-6)
  expect_gt(l2$p_value[l2$omitted == "a"], 0.1)
  expect_error(leave_one_out(bind_effects(eff("a", 0.1, 0.1))), "two studies")
  expect_error(leave_one_out(data.frame(a = 1)), "mm_effects")
})

# ================ meta-regression ================
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
  expect_error(meta_regress(eff, sd_, ~ nope), "not found in 'study_data'")
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

# ================ hierarchical FDR ================
pooled_tab <- function(feature, p, rank, part = "abundance", term = "x") {
  x <- data.frame(feature = feature, rank = rank, part = part, term = term, k = 3, estimate = 0.5, se = 0.1, ci_lo = 0, ci_hi = 1,
                  statistic = 1, p_value = p, tau2 = 0, I2 = 0, Q = 0, Q_p = 1, pi_lo = NA_real_, pi_hi = NA_real_,
                  max_weight = 0.4, top_study = "s1", p_adjust = NA_real_, stringsAsFactors = FALSE)
  class(x) <- c("mm_pooled", "data.frame"); x
}
feat <- c("P1", "P2", "C1a", "C1b", "C2a", "C2b", "G1", "G2", "G3")
par  <- c(P1 = NA, P2 = NA, C1a = "P1", C1b = "P1", C2a = "P2", C2b = "P2", G1 = "C1a", G2 = "C1a", G3 = "C2a")
pv   <- c(1e-4, 0.6, 0.001, 0.2, 0.0005, 0.3, 1e-5, 0.5, 1e-6)

test_that("hierarchical FDR tests only the children of rejected features, with the scaled level", {
  h <- adjust_hierarchical(pooled_tab(feat, pv, c("phylum", "phylum", rep("class", 4), rep("genus", 3))), par, q = 0.05)
  g <- stats::setNames(h$p_hier, h$feature)
  expect_equal(unname(h$level), c(1, 1, 2, 2, 2, 2, 3, 3, 3))
  # level 1: BH over the two phyla; one of two rejected, so level 2 is scaled by 1/2
  expect_equal(g[["P1"]], 2e-4); expect_equal(g[["P2"]], 0.6)
  expect_equal(g[["C1a"]], 0.002 / 0.5); expect_equal(g[["C1b"]], 0.2 / 0.5)
  expect_true(is.na(g[["C2a"]]) && is.na(g[["C2b"]]))                   # parent P2 not rejected: not tested
  # level 3 is scaled by (1/2) * (1/4): one of 2 phyla and one of 4 classes rejected
  expect_equal(g[["G1"]], 2e-5 / 0.125); expect_equal(g[["G2"]], min(0.5 / 0.125, 1))
  expect_true(is.na(g[["G3"]]))                                          # despite p = 1e-6: its ancestors were not rejected
  expect_equal(h$feature[h$hier_rejected], c("P1", "C1a", "G1"))
  expect_equal(attr(h, "settings")$hierarchical_q, 0.05)
})

test_that("nothing below an unrejected level is tested, and orphans become roots", {
  h <- adjust_hierarchical(pooled_tab(feat[1:4], c(0.5, 0.6, 1e-9, 1e-9), c("phylum", "phylum", "class", "class")), par, q = 0.05)
  expect_equal(sum(h$hier_rejected), 0); expect_true(all(is.na(h$p_hier[3:4])))
  # a class whose phylum is not in the table is treated as a root
  o <- adjust_hierarchical(pooled_tab(c("C1a", "C1b"), c(0.001, 0.2), c("class", "class")), par, q = 0.05)
  expect_equal(o$level, c(1, 1)); expect_equal(o$feature[o$hier_rejected], "C1a")
})

test_that("families are separate for each part and term", {
  a <- pooled_tab(c("P1", "P2"), c(1e-4, 0.6), "phylum", part = "abundance")
  b <- pooled_tab(c("P1", "P2"), c(0.7, 0.8), "phylum", part = "prevalence")
  h <- adjust_hierarchical(rbind(a, b), par)
  expect_equal(h$feature[h$hier_rejected & h$part == "abundance"], "P1"); expect_false(any(h$hier_rejected[h$part == "prevalence"]))
})

test_that("input is validated", {
  expect_error(adjust_hierarchical(data.frame(a = 1), par), "mm_pooled")
  expect_error(adjust_hierarchical(pooled_tab("P1", 0.1, "phylum"), unname(par)), "named")
  expect_error(adjust_hierarchical(pooled_tab("P1", 0.1, "phylum"), par, q = 2), "'q'")
})

test_that("feature_parents finds the nearest higher feature with the same lineage", {
  tax <- data.frame(kingdom = "Bacteria", phylum = c("Firmicutes", "Bacteroidota", "Firmicutes", "Firmicutes", "Bacteroidota", "Firmicutes"),
                    class = c(NA, NA, "Clostridia", "Clostridia", "Bacteroidia", "Bacilli"),
                    genus = c(NA, NA, NA, "Blautia", NA, "Lactobacillus"),
                    rank = c("phylum", "phylum", "class", "genus", "class", "genus"),
                    row.names = c("Fir", "Bac", "Clo", "Blau", "Bcl", "Lacto"), stringsAsFactors = FALSE)
  rel <- matrix(c(0.6, 0.4, 0.5, 0.3, 0.4, 0.2), 6, 2, dimnames = list(rownames(tax), c("S1", "S2")))
  s <- as_study(relabund = rel * 0.5, taxonomy = tax, samples = data.frame(g = 1:2, row.names = c("S1", "S2")), study = list(id = "h", feature_type = "taxon"))
  p <- feature_parents(s)
  expect_true(is.na(p[["Fir"]]) && is.na(p[["Bac"]]))
  expect_equal(p[["Clo"]], "Fir"); expect_equal(p[["Bcl"]], "Bac")
  expect_equal(p[["Blau"]], "Clo")                       # genus: class is the nearest present higher rank
  expect_equal(p[["Lacto"]], "Fir")                      # no Bacilli class feature: falls back to the phylum
})
