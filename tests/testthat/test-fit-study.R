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
  e <- fit_study(s, ~ group + age)
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
  e <- fit_study(make_bezi_study(), ~ group)
  expect_false("Rare" %in% e$feature)
  e2 <- fit_study(make_bezi_study(), ~ group, min_prevalence = 0.01)
  expect_true("Rare" %in% attr(e2, "fit_info")$feature)
})

test_that("a longitudinal model adds a subject random intercept", {
  skip_if_not_installed("gamlss")
  s <- make_bezi_study(subject = "person")
  e <- fit_study(s, ~ group)
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

test_that("fit_study reproduces taxa.compare() from 1.x on the Bangladesh data", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  utils::data("taxtab6", package = "metamicrobiomeR", envir = environment())
  tab <- as.data.frame(get("taxtab6", envir = environment()))
  phyla <- paste0("k__bacteria.p__", c("actinobacteria", "bacteroidetes", "firmicutes", "proteobacteria"))
  rel <- t(as.matrix(tab[, phyla])); colnames(rel) <- as.character(tab$x.sampleid)
  smp <- tab[, c("bf", "age.sample", "personid")]; rownames(smp) <- colnames(rel)
  s <- as_study(relabund = rel, samples = smp,
                study = list(id = "bangladesh", feature_type = "otu", subject = "personid"))
  new <- fit_study(s, ~ bf + age.sample)
  old <- suppressMessages(suppressWarnings(
    taxa.compare(taxtab = tab[, c("personid", "bf", "age.sample", phyla)], propmed.rel = "gamlss",
                 comvar = "bf", adjustvar = "age.sample", longitudinal = "yes")))
  for (term in c("bfNon_exclusiveBF", "bfNo_BF", "age.sample")) {
    n <- new[new$term == term, ]; n <- n[match(old$id, n$feature), ]
    expect_equal(n$estimate, old[[paste0("Estimate.", term)]], tolerance = 1e-8)
    expect_equal(n$se, old[[paste0("Std. Error.", term)]], tolerance = 1e-8)
    expect_equal(n$p_value, old[[paste0("Pr(>|t|).", term)]], tolerance = 1e-8)
  }
})
