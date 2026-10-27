eff <- function(study, scale = "logit_mean", feature = c("g1", "g2"), est = c(0.5, -0.2)) {
  new_effect_table(study = study, feature = feature, term = "groupy", estimate = est,
                   se = c(0.1, 0.2), part = "abundance", n = 20L, engine = "test",
                   effect_scale = scale)
}

test_that("new_effect_table recycles scalars and has the documented columns", {
  e <- eff("s1")
  expect_s3_class(e, "mm_effects")
  expect_s3_class(e, "data.frame")
  expect_equal(nrow(e), 2)
  expect_setequal(colnames(e), .mm_effect_cols)
  expect_true(all(e$converged))
})

test_that("check_effects rejects invalid tables", {
  expect_error(eff("s1", est = c(Inf, 1)), "finite")
  expect_error(new_effect_table("s", "f", "t", 1, -0.1), "negative")
  expect_error(new_effect_table("s", "f", "t", 1, 0.1, part = "weird"), "part")
  expect_error(new_effect_table(c("s", "s"), c("f", "f"), "t", c(1, 2), c(0.1, 0.1)), "Duplicate")
  expect_warning(new_effect_table("s", "f", "t", 1, 0), "zero")
})

test_that("bind_effects stacks studies and refuses mixed effect scales", {
  b <- bind_effects(eff("s1"), eff("s2"))
  expect_equal(nrow(b), 4)
  expect_equal(sort(unique(b$study)), c("s1", "s2"))
  expect_equal(nrow(bind_effects(list(eff("s1"), eff("s2")))), 4)
  expect_error(bind_effects(eff("s1"), eff("s2", scale = "smd")), "effect scales")
  expect_equal(nrow(bind_effects(eff("s1"), eff("s2", scale = "smd"), allow_mixed_scale = TRUE)), 4)
  expect_error(bind_effects(eff("s1"), eff("s1")), "Duplicate")
  expect_error(bind_effects(data.frame(a = 1)), "mm_effects")
})
