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
  p <- pool_effects(bound, method = "DL", hk = FALSE)
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
  a <- pool_effects(bound, method = "DL", hk = FALSE)
  b <- pool_effects(bound, method = "DL", hk = TRUE)
  expect_equal(a$estimate, b$estimate, tolerance = 1e-8)
  expect_gt(b$ci_hi - b$ci_lo, a$ci_hi - a$ci_lo)
  expect_false(isTRUE(all.equal(a$p_value, b$p_value)))
  expect_false(is.na(b$pi_lo))              # prediction interval with 4 studies
  two <- pool_effects(bind_effects(eff("s1", 0.5, 0.2), eff("s2", 0.4, 0.2)), method = "DL")
  expect_true(is.na(two$pi_lo))             # not reported with fewer than 3 studies
})

test_that("features need enough studies, and a missing feature is not pooled", {
  b <- bind_effects(eff("s1", 0.5, 0.2, "g1"), eff("s2", 0.4, 0.2, "g1"), eff("s3", 0.6, 0.2, "g1"),
                    eff("s1", 0.1, 0.2, "g2"),                       # only one study
                    eff("s2", 0.3, 0.2, "g3"), eff("s3", 0.2, 0.2, "g3"))
  p <- pool_effects(b, method = "DL", hk = FALSE)
  expect_setequal(p$feature, c("g1", "g3"))
  expect_equal(p$k[p$feature == "g3"], 2)
  p2 <- pool_effects(b, method = "DL", hk = FALSE, min_prop = 1)       # must be in all 3 studies
  expect_equal(p2$feature, "g1")
})

test_that("p-values are adjusted within the requested family", {
  mk <- function(feature, rank, est) bind_effects(eff("s1", est, 0.1, feature, rank = rank), eff("s2", est, 0.1, feature, rank = rank))
  b <- bind_effects(mk("a", "genus", 0.5), mk("b", "genus", 0.1), mk("c", "family", 0.5), mk("d", "family", 0.1))
  by_rank <- pool_effects(b, method = "DL", hk = FALSE)                                   # default: within rank
  mixed <- pool_effects(b, method = "DL", hk = FALSE, adjust_by = c("part", "term"))
  expect_equal(by_rank$p_adjust[by_rank$feature == "b"], stats::p.adjust(by_rank$p_value[by_rank$rank == "genus"], "fdr")[by_rank$feature[by_rank$rank == "genus"] == "b"])
  expect_equal(sort(mixed$p_adjust), sort(stats::p.adjust(mixed$p_value, "fdr")))
  expect_equal(pool_effects(b, method = "DL", hk = FALSE, p_adjust = "none")$p_adjust,
               pool_effects(b, method = "DL", hk = FALSE, p_adjust = "none")$p_value)
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
  expect_warning(p <- pool_effects(b, method = "DL"), "dropped")
  expect_equal(p$k, 2)
  expect_error(pool_effects(bind_effects(eff("s1", 0.5, 0.2)), min_studies = 2), "No feature")
})
