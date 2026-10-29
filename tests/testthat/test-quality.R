mkq <- function(study, est, se, feature = paste0("f", seq_along(est)), converged = TRUE) {
  new_effect_table(study = study, feature = feature, term = "x", estimate = est, se = se, rank = "genus", part = "abundance",
                   n = 20L, engine = "t", effect_scale = "logit_mean", converged = converged)
}
set.seed(2)
base_se <- function(n) runif(n, 0.15, 0.45)
A <- mkq("A", rnorm(10, 0.3, 0.2), base_se(10)); B <- mkq("B", rnorm(10, 0.3, 0.2), base_se(10))
C <- mkq("C", c(rnorm(7, 0.3, 0.2), 0.43, -10, 3.7), c(base_se(7), 0.005, 0.003, 0.009))        # three degenerate fits, like USA(NC)
eff <- bind_effects(A, B, C)

test_that("check_effects_quality flags implausibly small standard errors of one study", {
  q <- check_effects_quality(eff)
  expect_s3_class(q, "mm_quality")
  r <- q$rows
  expect_equal(sum(r$flag_se), 3); expect_true(all(r$study[r$flag_se] == "C")); expect_equal(sum(r$flag_se[r$study != "C"]), 0)
  expect_equal(q$studies$flagged_se[q$studies$study == "C"], 3); expect_lt(q$studies$se_min[q$studies$study == "C"], 0.01)
  expect_gt(q$studies$reference_se_q[q$studies$study == "C"], 0.1)
  expect_output(print(q), "study C")
  expect_equal(sum(check_effects_quality(eff, se_ratio = 0.001)$rows$flag_se), 0)                  # a stricter ratio flags nothing
  big <- check_effects_quality(eff, max_abs_estimate = 5)$rows; expect_equal(sum(big$flag_estimate), 1); expect_equal(abs(big$estimate[big$flag_estimate]), 10)
  nc <- check_effects_quality(bind_effects(A, mkq("B", rnorm(10), base_se(10), converged = c(rep(TRUE, 8), FALSE, FALSE))))
  expect_equal(sum(nc$rows$flag_converged), 2)
  one <- check_effects_quality(A); expect_equal(sum(one$rows$flag_se), 0)                          # no reference without other studies
  expect_error(check_effects_quality(data.frame(a = 1)), "mm_effects")
})

test_that("filter_effects removes the flagged estimates, records why, and changes the pooled result", {
  expect_message(f <- filter_effects(eff), "3 of 30 estimates dropped")
  expect_s3_class(f, "mm_effects"); expect_equal(nrow(f), 27)
  dr <- attr(f, "dropped"); expect_equal(nrow(dr), 3); expect_true(all(grepl("implausible standard error", dr$reason)))
  # one feature whose pooled result depends on the degenerate study
  dom <- bind_effects(mkq("A", 0.01, 0.43, "g"), mkq("B", 0.55, 0.35, "g"), mkq("C", 0.43, 0.005, "g"), mkq("A", 0.2, 0.3, "h"), mkq("B", 0.3, 0.3, "h"), mkq("C", 0.25, 0.3, "h"))
  before <- pool_effects(dom, method = "DL", test = "z"); after <- pool_effects(suppressMessages(filter_effects(dom)), method = "DL", test = "z")
  expect_lt(before$p_value[before$feature == "g"], 1e-6)
  expect_equal(after$k[after$feature == "g"], 2); expect_gt(after$p_value[after$feature == "g"], 0.1)
  expect_equal(nrow(suppressMessages(filter_effects(eff, drop_implausible_se = FALSE))), 30)
  bad <- bind_effects(A, mkq("B", rnorm(10, 0.3, 0.2), base_se(10), converged = c(rep(TRUE, 9), FALSE)))
  expect_equal(nrow(suppressMessages(filter_effects(bad))), 19); expect_equal(nrow(filter_effects(bad, drop_nonconverged = FALSE)), 20)
  expect_silent(filter_effects(A))
})

test_that("pool_sensitivity shows settings that disagree", {
  near <- bind_effects(mkq("s1", 0.050, 0.5, "near_identical"), mkq("s2", 0.052, 0.5, "near_identical"))
  real <- bind_effects(mkq(paste0("s", 1:5), c(0.50, 0.60, 0.55, 0.45, 0.50), rep(0.1, 5), "real"))
  two <- bind_effects(near, real)
  s <- pool_sensitivity(two, significance = "p_value", min_prop = 0)     # near_identical is in 2 of 5 studies
  expect_s3_class(s, "mm_poolcompare")
  expect_setequal(names(s$n_significant), c("DL_z", "REML_z", "REML_knha", "REML_adhoc"))
  w <- s$wide
  ni <- w[w$feature == "near_identical", ]
  expect_true(ni$sig_REML_knha); expect_false(ni$sig_DL_z); expect_false(ni$sig_REML_adhoc)      # plain Hartung-Knapp is the odd one out
  expect_equal(ni$n_settings_significant, 1); expect_false(ni$robust)
  expect_true(w$robust[w$feature == "real"]); expect_equal(s$robust, "real")
  expect_equal(nrow(s$long), 8)
  expect_output(print(s), "significant under all 4 settings: 1 of 2")
  # with only two studies the Hartung-Knapp test has one degree of freedom: even a clear effect is not significant
  clear2 <- bind_effects(mkq("a", 0.5, 0.1, "g"), mkq("b", 0.6, 0.1, "g"))
  c2 <- pool_sensitivity(clear2, significance = "p_value")$wide
  expect_true(c2$sig_DL_z); expect_false(c2$sig_REML_adhoc)
  custom <- pool_sensitivity(two, settings = list(a = list(method = "DL", test = "z"), b = list(method = "PM", test = "z")), significance = "p_value", min_prop = 0)
  expect_equal(names(custom$n_significant), c("a", "b"))
  expect_error(pool_sensitivity(two, settings = list(list(method = "DL"))), "named list")
})

test_that("combine_pvalues matches the hand calculation and reports direction and sign test", {
  e <- bind_effects(new_effect_table("s1", "f", "x", 0.5, 0.2, p_value = 0.01, n = 30L, engine = "a", effect_scale = "s"),
                    new_effect_table("s2", "f", "x", 0.3, 0.2, p_value = 0.04, n = 50L, engine = "b", effect_scale = "t"),
                    new_effect_table("s3", "f", "x", -0.1, 0.2, p_value = 0.20, n = 40L, engine = "c", effect_scale = "u"), allow_mixed_scale = TRUE, allow_mixed_engine = TRUE)
  cp <- combine_pvalues(e)                                                                         # different engines and scales are fine
  w <- sqrt(c(30, 50, 40)); z <- c(1, 1, -1) * stats::qnorm(c(0.01, 0.04, 0.20) / 2, lower.tail = FALSE)
  Z <- sum(w * z) / sqrt(sum(w^2))
  expect_s3_class(cp, "mm_combined")
  expect_equal(cp$statistic, Z, tolerance = 1e-10); expect_equal(cp$p_value, 2 * stats::pnorm(-abs(Z)), tolerance = 1e-10)
  expect_equal(c(cp$k, cp$n_pos, cp$n_neg), c(3, 2, 1)); expect_equal(cp$direction, "mixed"); expect_equal(cp$sign_test_p, 1)
  unw <- combine_pvalues(e, weights = "none"); expect_equal(unw$statistic, sum(z) / sqrt(3), tolerance = 1e-10)
  fi <- combine_pvalues(e, method = "fisher")
  pp <- c(0.005, 0.02, 1 - 0.10); pn <- 1 - pp
  expect_equal(fi$p_value, min(1, 2 * min(stats::pchisq(-2 * sum(log(pp)), 6, lower.tail = FALSE), stats::pchisq(-2 * sum(log(pn)), 6, lower.tail = FALSE))), tolerance = 1e-10)
  cons <- bind_effects(mkq("a", 0.4, 0.1, "g"), mkq("b", 0.5, 0.1, "g"), mkq("c", 0.6, 0.1, "g"))
  cc <- combine_pvalues(cons); expect_equal(cc$direction, "up"); expect_lt(cc$p_value, 1e-6); expect_equal(cc$sign_test_p, 0.25)    # 3 of 3 positive: p = 2 * 0.5^3
})

test_that("combine_pvalues computes missing p-values, adjusts within the family and validates", {
  e <- bind_effects(mkq("a", c(0.5, 0.0), c(0.1, 0.1), c("f", "g")), mkq("b", c(0.4, 0.1), c(0.1, 0.1), c("f", "g")))
  cp <- combine_pvalues(e)                                                       # p_value is NA: derived from estimate and se
  expect_equal(nrow(cp), 2); expect_true(all(!is.na(cp$p_value))); expect_equal(cp$p_adjust, stats::p.adjust(cp$p_value, "fdr"))
  expect_equal(nrow(combine_pvalues(e, min_studies = 2, term = "x")), 2)
  expect_error(combine_pvalues(e, min_studies = 3), "at least 3 studies"); expect_error(combine_pvalues(data.frame(a = 1)), "mm_effects")
})
