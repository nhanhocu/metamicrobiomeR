make_count_study <- function(n = 100, seed = 21) {
  set.seed(seed)
  grp <- factor(rep(c("a", "b"), length.out = n)); age <- runif(n, 0, 6)
  N <- round(exp(stats::rnorm(n, log(20000), 0.6)))            # depth varies several-fold
  sim <- function(mu, phi = 10, zi = 0) {
    p <- stats::rbeta(n, mu * phi, (1 - mu) * phi); y <- stats::rbinom(n, N, p)
    y[stats::rbinom(n, 1, zi) == 1] <- 0; y
  }
  counts <- rbind(F1 = sim(stats::plogis(-4 + 0.9 * (grp == "b"))),                       # abundance effect 0.9
                  F2 = sim(rep(stats::plogis(-4), n)),                                   # null
                  F3 = sim(rep(stats::plogis(-4), n), zi = ifelse(grp == "b", 0.5, 0.05)),  # prevalence effect
                  F4 = sim(rep(stats::plogis(-1.5), n)))                                 # abundant: no zeros
  colnames(counts) <- paste0("S", seq_len(n))
  smp <- data.frame(group = grp, age = age, person = rep(seq_len(n / 2), each = 2), row.names = colnames(counts))
  as_study(counts = counts, samples = smp, study = list(id = "cnt", feature_type = "asv", subject = "person"))
}

test_that("zibb_glmmtmb recovers an abundance effect and uses the sequencing depth", {
  skip_if_not_installed("glmmTMB")
  s <- make_count_study()
  e <- fit_study(s, ~ group + age, engine = "zibb_glmmtmb", longitudinal = FALSE)
  expect_s3_class(e, "mm_effects")
  expect_true(all(e$engine == "zibb_glmmtmb" & e$effect_scale == "logit_mean" & e$part == "abundance"))
  g <- e[e$feature == "F1" & e$term == "groupb", ]
  expect_gt(g$estimate, 0.4); expect_lt(g$estimate, 1.5); expect_lt(g$p_value, 0.05)
  nul <- e[e$feature == "F2" & e$term == "groupb", ]
  expect_lt(abs(nul$estimate), 0.8)
  expect_true(all(attr(e, "fit_info")$converged))
  # an all-zero-depth sample is dropped, not fitted
  cnt <- mm_counts(s); cnt[, 5] <- 0L
  s0 <- suppressWarnings(as_study(counts = cnt, samples = mm_samples(s), study = list(id = "z", feature_type = "asv")))
  e0 <- fit_study(s0, ~ group, engine = "zibb_glmmtmb", longitudinal = FALSE)
  expect_equal(unique(e0$n), 99L)
})

test_that("zibb_glmmtmb needs counts and supports a subject random intercept", {
  skip_if_not_installed("glmmTMB")
  s <- make_count_study()
  rel <- as_study(relabund = mm_relabund(s), samples = mm_samples(s), study = list(id = "r", feature_type = "taxon"))
  expect_error(fit_study(rel, ~ group, engine = "zibb_glmmtmb"), "needs counts")
  e <- fit_study(s, ~ group, engine = "zibb_glmmtmb")             # subject is set: mixed model
  expect_true(all(e$term == "groupb"))
  expect_gt(e$estimate[e$feature == "F1"], 0.3)
})

test_that("the prevalence part is opt-in, labelled, and unidentifiable estimates are flagged", {
  skip_if_not_installed("glmmTMB")
  s <- make_count_study()
  base <- fit_study(s, ~ group, engine = "zibb_glmmtmb", longitudinal = FALSE)
  expect_false("prevalence" %in% base$part)
  e <- fit_study(s, ~ group, engine = "zibb_glmmtmb", longitudinal = FALSE, prevalence = ~ group)
  prev <- e[e$part == "prevalence", ]
  expect_gt(nrow(prev), 0)
  expect_true(all(prev$effect_scale == "logit_zero_prob"))
  expect_false("F4" %in% prev$feature)                            # no zeros: zero-inflation part dropped
  f3 <- prev[prev$feature == "F3" & prev$term == "groupb", ]
  expect_gt(f3$estimate, 0)                                        # more structural zeros in group b
  expect_equal(f3$converged, is.finite(f3$se) && f3$se <= 10 && attr(e, "fit_info")$converged[attr(e, "fit_info")$feature == "F3"])
  plain <- fit_study(s, ~ group, engine = "zibb_glmmtmb", longitudinal = FALSE, zero_inflation = FALSE)
  expect_true(all(plain$part == "abundance"))
  expect_error(fit_study(s, ~ group, engine = "zibb_glmmtmb", prevalence = "group"), "one-sided")
})

test_that("lm_clr works on the centred log-ratio of the analysed features", {
  s <- make_count_study()
  clr <- metamicrobiomeR:::.clr_matrix(s, seq_len(nrow(mm_counts(s))), NULL)
  expect_equal(unname(colSums(clr)), rep(0, ncol(clr)), tolerance = 1e-8)         # centred
  e <- fit_study(s, ~ group, engine = "lm_clr", longitudinal = FALSE, min_prevalence = 0)
  expect_true(all(e$engine == "lm_clr" & e$effect_scale == "clr_mean"))
  # F1 rises with group b relative to the other features
  expect_gt(e$estimate[e$feature == "F1" & e$term == "groupb"], 0.3)
  # the composition is the full set of features at the rank, whatever passes the filter
  e2 <- fit_study(s, ~ group, engine = "lm_clr", longitudinal = FALSE, min_prevalence = 0.99)
  expect_equal(e2$estimate[e2$feature == "F4"], e$estimate[e$feature == "F4"], tolerance = 1e-10)
})

test_that("estimates from different engines are not combined by bind_effects", {
  a <- new_effect_table("s1", "f", "x", 0.5, 0.1, engine = "zibb_glmmtmb", effect_scale = "logit_mean")
  b <- new_effect_table("s2", "f", "x", 0.4, 0.1, engine = "bezi_gamlss", effect_scale = "logit_mean")
  expect_error(bind_effects(a, b), "different engines")
  expect_equal(nrow(bind_effects(a, b, allow_mixed_engine = TRUE)), 2)
})
