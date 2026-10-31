# Tests for the meta-analyse family; merged from several smaller files (see git history).

# ================ meta_analyse ================
sims <- function(n = 3, ...) simulate_studies(n_studies = n, n_samples = 30, n_features = 8, n_diff = 3, effect = 1.5, tau = 0.2, seed = 4, ...)$studies
run <- function(s, ...) suppressMessages(meta_analyse(s, ~ group, term = "groupb", engine = "lm_asin", fit_args = list(longitudinal = FALSE, min_prevalence = 0), ...))

test_that("meta_analyse equals the step-by-step workflow and returns one object", {
  s <- sims()
  res <- run(s)
  expect_s3_class(res, "mm_meta")
  manual <- pool_effects(bind_effects(lapply(s, fit_study, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0)), term = "groupb")
  expect_equal(as.data.frame(res$pooled), as.data.frame(manual))
  expect_equal(as.data.frame(res), as.data.frame(res$pooled))
  expect_s3_class(res$effects, "mm_effects"); expect_equal(res$studies$study, c("sim1", "sim2", "sim3"))
  expect_equal(res$settings$engine, "lm_asin"); expect_equal(res$settings$formula, "~group")
  out <- testthat::capture_output(print(res))
  expect_match(out, "<mm_meta> 3 studies \\(sim1, sim2, sim3\\) \\| engine lm_asin \\| model ~group")
  expect_match(out, "8 features pooled \\(REML, adhoc test\\) \\| adjusted p <= 0.05:"); expect_match(out, "checks: 0 warning\\(s\\)")
  expect_match(out, "plot\\(x\\) draws a forest plot")
})

test_that("the summary explains an empty result instead of showing a bare zero", {
  few <- simulate_studies(n_studies = 3, n_samples = 14, n_features = 6, n_diff = 1, effect = 0.2, tau = 0.5, seed = 2)$studies
  res <- suppressMessages(meta_analyse(few, ~ group, term = "groupb", engine = "lm_asin", fit_args = list(longitudinal = FALSE, min_prevalence = 0)))
  expect_match(testthat::capture_output(print(res)), "No feature reaches adjusted p <= 0.05; the table is ordered by p-value. With 3 studies there is little power")
  strong <- suppressMessages(meta_analyse(simulate_studies(n_studies = 5, n_samples = 80, n_features = 6, n_diff = 3, effect = 3, tau = 0.05, depth_sd = 0.2, seed = 3)$studies,
                                          ~ group, term = "groupb", engine = "lm_asin", test = "z", fit_args = list(longitudinal = FALSE, min_prevalence = 0)))
  expect_false(grepl("No feature reaches", testthat::capture_output(print(strong))))
})

test_that("plot and as.data.frame work from the result", {
  res <- run(sims())
  f <- plot(res); expect_s3_class(f, "ggplot")
  expect_length(unique(ggplot2::ggplot_build(f)$layout$layout$PANEL), 3)                          # the three strongest features
  expect_length(unique(ggplot2::ggplot_build(plot(res, top = 2))$layout$layout$PANEL), 2)
  one <- plot(res, feature = "F1"); expect_length(unique(ggplot2::ggplot_build(one)$layout$layout$PANEL), 1)
  expect_s3_class(plot(res, "heatmap"), "ggplot"); expect_s3_class(plot(res, "volcano"), "ggplot")
  with_x <- lapply(sims(), function(st) { set.seed(1); st$samples$x <- rnorm(nrow(st$samples)); st })
  both <- suppressMessages(meta_analyse(with_x, ~ group + x, engine = "lm_asin", fit_args = list(longitudinal = FALSE, min_prevalence = 0)))
  expect_error(plot(both), "Several terms were pooled"); expect_s3_class(plot(both, term = "groupb"), "ggplot")
})

test_that("the engine is chosen automatically and the choice is reported", {
  s <- sims()
  skip_if_not_installed("glmmTMB")
  m1 <- testthat::capture_messages(meta_analyse(s, ~ group, term = "groupb", fit_args = list(longitudinal = FALSE, min_prevalence = 0.5)))
  expect_true(any(grepl("Using the 'zibb_glmmtmb' engine \\(every study has counts; about a second per feature and study", m1)))
  rel <- lapply(s, function(x) as_study(relabund = mm_relabund(x), samples = mm_samples(x), study = list(id = x$study$id, feature_type = "taxon")))
  m2 <- testthat::capture_messages(meta_analyse(rel, ~ group, term = "groupb", fit_args = list(longitudinal = FALSE, min_prevalence = 0)))
  expect_true(any(grepl("Using the 'lm_asin' engine \\(not every study has counts\\)", m2)))
})

# a study whose fits are degenerate: within each group the abundances are (almost) constant, so the standard errors collapse
degenerate <- function(base) {
  set.seed(9); rel <- mm_relabund(base); g <- mm_samples(base)$group
  for (j in seq_len(nrow(rel))) for (lv in levels(g)) rel[j, g == lv] <- mean(rel[j, g == lv]) + rnorm(sum(g == lv), 0, 1e-9)
  rel <- pmax(rel, 1e-6); rel <- sweep(rel, 2, colSums(rel), "/") * 0.999
  as_study(relabund = rel, samples = mm_samples(base), study = list(id = "bad", feature_type = "taxon"))
}

test_that("implausible estimates are reported by default and can be filtered out", {
  s <- c(sims(2), list(degenerate(sims(1)[[1]])))
  expect_warning(res <- suppressMessages(meta_analyse(s, ~ group, term = "groupb", engine = "lm_asin", fit_args = list(longitudinal = FALSE, min_prevalence = 0))),
                 "estimates look doubtful \\(bad: \\d+\\).*use quality = \"filter\" to drop them")
  expect_gt(sum(res$quality$rows$flagged), 0); expect_equal(nrow(res$effects), 24)                                # kept
  expect_message(f <- meta_analyse(s, ~ group, term = "groupb", engine = "lm_asin", quality = "filter", fit_args = list(longitudinal = FALSE, min_prevalence = 0)),
                 "Dropped \\d+ of 24 estimates before pooling")
  expect_lt(nrow(f$effects), 24)
  expect_silent(suppressMessages(meta_analyse(s, ~ group, term = "groupb", engine = "lm_asin", quality = "off", fit_args = list(longitudinal = FALSE, min_prevalence = 0))))
  expect_null(suppressWarnings(suppressMessages(meta_analyse(s, ~ group, term = "groupb", engine = "lm_asin", quality = "off", fit_args = list(longitudinal = FALSE, min_prevalence = 0))))$quality)
})

test_that("warnings from a study are passed on with its name, not lost", {
  s <- sims(); one <- s[[2]]; one$samples$group <- factor(rep("a", 30)); s[[2]] <- one                  # every fit in sim2 fails
  expect_warning(suppressMessages(meta_analyse(s, ~ group, term = "groupb", engine = "lm_asin", fit_args = list(longitudinal = FALSE, min_prevalence = 0), quality = "off")),
                 "^\\[sim2\\] \\d+ of 8 fits failed")
})

test_that("problems are reported with the name of the study they come from", {
  s <- sims(); na <- s[[2]]; na$samples$group[1:5] <- NA; s[[2]] <- na
  m <- testthat::capture_messages(meta_analyse(s, ~ group, term = "groupb", engine = "lm_asin", fit_args = list(longitudinal = FALSE, min_prevalence = 0)))
  expect_true(any(grepl("^\\[sim2\\] 5 of 30 samples have missing values", m)))
  bad <- sims(); bad[[3]]$samples$group <- NULL
  expect_error(meta_analyse(bad, ~ group, engine = "lm_asin"), "Study 'sim3': Variable not found in the sample data: group")
  small <- c(sims(2), list(simulate_study(n_samples = 8, n_features = 8, n_diff = 3, id = "tiny", seed = 3)$study))
  expect_warning(suppressMessages(meta_analyse(small, ~ group, term = "groupb", engine = "lm_asin", fit_args = list(longitudinal = FALSE, min_prevalence = 0))),
                 "The studies may not be safe to pool: small studies \\(fewer than 10 samples\\) \\(tiny \\(n = 8\\)\\)")
  expect_silent(suppressMessages(suppressWarnings(meta_analyse(small, ~ group, term = "groupb", engine = "lm_asin", check = FALSE, quality = "off", fit_args = list(longitudinal = FALSE, min_prevalence = 0)))))
})

test_that("a rank agglomerates ASV-level studies and matches taxa across databases", {
  mk <- function(id, genera, phyla, seed) {
    set.seed(seed); n <- 20; nf <- length(genera)
    cnt <- matrix(rpois(nf * n, 40), nf, n, dimnames = list(paste0("f", seq_len(nf)), paste0("S", seq_len(n))))
    cnt[1, rep(c(FALSE, TRUE), 10)] <- cnt[1, rep(c(FALSE, TRUE), 10)] * 4L                              # one differential taxon
    as_study(counts = cnt, taxonomy = data.frame(phylum = phyla, genus = genera, row.names = rownames(cnt), stringsAsFactors = FALSE),
             samples = data.frame(group = factor(rep(c("a", "b"), 10)), row.names = colnames(cnt)), study = list(id = id, feature_type = "asv"))
  }
  a <- mk("A", c("Blautia", "Roseburia", "[Ruminococcus]", "Bacteroides"), c("F", "F", "F", "B"), 1)
  b <- mk("B", c("Blautia", "Ruminococcus", "Bacteroides", "Prevotella"), c("F_A", "F_A", "B", "B"), 2)
  c3 <- mk("C", c("Blautia", "Ruminococcus", "Bacteroides", "Alistipes"), c("F", "F", "B", "B"), 3)
  res <- suppressMessages(meta_analyse(list(a, b, c3), ~ group, term = "groupb", rank = "genus", engine = "lm_asin", fit_args = list(longitudinal = FALSE, min_prevalence = 0), quality = "off"))
  expect_true(all(c("Blautia", "Ruminococcus", "Bacteroides") %in% res$pooled$feature))                       # [Ruminococcus] matched with Ruminococcus
  expect_false("Prevotella" %in% res$pooled$feature)                                                          # in one study only
  expect_equal(unique(res$pooled$rank), "genus"); expect_equal(res$pooled$k[res$pooled$feature == "Blautia"], 3)
  expect_equal(res$settings$rank, "genus")
})

test_that("meta_analyse validates its input", {
  s <- sims()
  expect_error(meta_analyse(s[1], ~ group), "needs at least two studies; got 1. Use fit_study\\(\\)")
  expect_error(meta_analyse(list(s[[1]], s[[1]]), ~ group), "Study IDs must be unique; 'sim1' is used twice")
  expect_error(meta_analyse(s, group ~ 1), "'formula' must be one-sided")
  expect_error(meta_analyse(list(1, 2), ~ group), "Supply mm_study objects")
  expect_error(meta_analyse(s, ~ group, engine = "nonsense"))
})

test_that("summary describes between-study variation, replicability and what can explain it", {
  mk <- function(study, est, se, feature) new_effect_table(study, feature, "x", est, se, rank = "genus", part = "abundance", n = 30L, engine = "t", effect_scale = "s")
  eff <- bind_effects(mk(paste0("s", 1:5), c(0.1, 0.9, 0.5, 0.4, 0.6), rep(0.05, 5), "het"),      # present on average, studies disagree
                      mk(paste0("s", 1:5), c(0.50, 0.52, 0.49, 0.51, 0.50), rep(0.05, 5), "same"),  # studies agree
                      mk(paste0("s", 1:5), c(0.00, 0.02, -0.01, 0.01, 0.00), rep(0.05, 5), "null"))
  pooled <- pool_effects(eff)
  fake <- structure(list(pooled = pooled, effects = eff, studies = data.frame(study = paste0("s", 1:5)), checks = NULL, quality = NULL,
                         settings = list(formula = "~x", engine = "t", rank = NULL, term = NULL, method = "REML", test = "adhoc", quality = "off")), class = "mm_meta")
  s <- summary(fake)
  expect_s3_class(s, "summary.mm_meta")
  p <- as.data.frame(pooled)
  expect_equal(s$I2[2], stats::median(p$I2)); expect_equal(s$n_features, 3)
  expect_equal(s$n_excl0, 2); expect_equal(s$n_not_replicable, 1)                         # "het" and "same" excluded zero; only "het" has a prediction interval that includes it
  expect_equal(as.character(s$top$feature[1]), "het")
  out <- testthat::capture_output(print(s))
  expect_match(out, "I2 \\(share of variation beyond sampling error\\): median")
  expect_match(out, "for 1 the prediction interval still includes zero: present on average, but a new study could show no effect")
  expect_match(out, "Most heterogeneous:"); expect_match(out, "meta_regress\\(\\)"); expect_match(out, "~ group \\+ batch")
})

test_that("a high-heterogeneity hint appears in the printed summary", {
  s <- simulate_studies(n_studies = 5, n_samples = 60, n_features = 10, n_diff = 3, effect = 1.5, tau = 1.5, seed = 8)$studies
  res <- suppressMessages(meta_analyse(s, ~ group, term = "groupb", engine = "lm_asin", quality = "off", check = FALSE, fit_args = list(longitudinal = FALSE, min_prevalence = 0)))
  hi <- stats::median(as.data.frame(res)$I2)
  out <- testthat::capture_output(print(res))
  if (hi > 50) expect_match(out, "Heterogeneity is high \\(median I2 \\d+%\\): the studies disagree. summary\\(x\\) shows how much") else expect_false(grepl("Heterogeneity is high", out))
  quiet <- suppressMessages(meta_analyse(simulate_studies(n_studies = 4, n_samples = 40, n_features = 6, n_diff = 2, effect = 1.5, tau = 0, seed = 3)$studies, ~ group, term = "groupb",
                                         engine = "lm_asin", quality = "off", check = FALSE, fit_args = list(longitudinal = FALSE, min_prevalence = 0)))
  expect_false(grepl("Heterogeneity is high", testthat::capture_output(print(quiet))))
  expect_match(testthat::capture_output(print(quiet)), "summary\\(x\\) describes the between-study variation")
})

# ================ diversity and report ================
sims <- function(n = 3) simulate_studies(n_studies = n, n_samples = 30, n_features = 8, n_diff = 3, effect = 1.5, tau = 0.2, seed = 4)$studies
quiet <- function(...) suppressMessages(meta_analyse(..., engine = "lm_asin", fit_args = list(longitudinal = FALSE, min_prevalence = 0)))

test_that("meta_analyse says which study it is on", {
  m <- testthat::capture_messages(meta_analyse(sims(), ~ group, term = "groupb", engine = "lm_asin", fit_args = list(longitudinal = FALSE, min_prevalence = 0)))
  expect_true(all(c("Study 1 of 3: sim1", "Study 2 of 3: sim2", "Study 3 of 3: sim3") %in% trimws(m)))
  old <- options(metamicrobiomeR.quiet = TRUE); on.exit(options(old))
  expect_silent(meta_analyse(sims(), ~ group, term = "groupb", engine = "lm_asin", quality = "off", check = FALSE, fit_args = list(longitudinal = FALSE, min_prevalence = 0)))
})

test_that("alpha and beta diversity are available from meta_analyse and equal the step-by-step results", {
  skip_if_not_installed("vegan")
  s <- sims()
  res <- quiet(s, ~ group, term = "groupb", alpha = TRUE, beta = TRUE, beta_args = list(permutations = 99, seed = 1))
  alpha_manual <- pool_effects(bind_effects(lapply(s, fit_alpha, ~ group, longitudinal = FALSE)), term = "groupb", min_prop = 0)
  expect_equal(as.data.frame(res$alpha), as.data.frame(alpha_manual)); expect_setequal(res$alpha$feature, c("shannon", "simpson", "observed"))
  beta_manual <- pool_beta(lapply(s, fit_beta, ~ group, permutations = 99, seed = 1))
  expect_equal(res$beta, beta_manual); expect_s3_class(res$alpha_effects, "mm_effects")
  out <- testthat::capture_output(print(res))
  expect_match(out, "Alpha diversity \\(pooled standardised difference\\): .*shannon"); expect_match(out, "Beta diversity \\(PERMANOVA, combined across studies\\): group R2")
  plain <- quiet(s, ~ group, term = "groupb"); expect_null(plain$alpha); expect_null(plain$beta)
  expect_false(grepl("Alpha diversity|Beta diversity", testthat::capture_output(print(plain))))
  bad <- sims(); bad[[2]]$study$subject <- "person"; bad[[2]]$samples$person <- rep(1:15, each = 2)
  expect_error(quiet(bad, ~ group, term = "groupb", beta = TRUE, beta_args = list(permutations = 9)), "Study 'sim2' \\(beta diversity\\)")
})

test_that("write_report writes a readable report with figures", {
  skip_if_not_installed("vegan")
  res <- quiet(sims(), ~ group, term = "groupb", alpha = TRUE, beta = TRUE, beta_args = list(permutations = 49, seed = 1))
  f <- file.path(tempfile(fileext = ".dir"), "report.md"); dir.create(dirname(f))
  expect_identical(write_report(res, f, top = 3), f)
  r <- readLines(f)
  expect_true(all(c("# Meta-analysis report", "## Analysis", "## Studies", "## Checks before pooling", "## Between-study variation", "## Alpha diversity", "## Beta diversity", "## Figures", "## Software") %in% r))
  expect_true(any(grepl("^## Pooled results \\(top 3 of 8 features", r))); expect_true(any(grepl("^- Model: `~group`; engine: `lm_asin`", r)))
  expect_true(any(grepl("^\\| sim1 \\|", r)))                                                      # the studies table
  p <- as.data.frame(res); top <- p$feature[which.min(p$p_value)]; expect_true(any(grepl(paste0("^\\| ", top, " \\|"), r)))
  expect_true(any(grepl("!\\[forest\\]\\(report_files/forest.png\\)", r)))
  expect_true(all(file.exists(file.path(dirname(f), "report_files", c("forest.png", "volcano.png", "heatmap.png")))))
  expect_true(any(grepl("Significant under all settings", r))); expect_true(any(grepl("metafor: ", r)))
  g <- file.path(dirname(f), "nofig.md"); write_report(res, g, figures = FALSE); expect_false(any(grepl("## Figures", readLines(g))))
  expect_error(write_report(1), "'x' must be an mm_meta object")
})

test_that("write_report copes with several pooled terms", {
  with_x <- lapply(sims(), function(st) { set.seed(1); st$samples$x <- rnorm(nrow(st$samples)); st })
  res <- quiet(with_x, ~ group + x)
  f <- tempfile(fileext = ".md"); expect_silent(write_report(res, f))
  r <- readLines(f)
  expect_true(any(grepl("Not drawn: forest, heatmap\\. This happens when several terms were pooled", r)))     # the volcano plot does not need a single term
  expect_true(any(grepl("!\\[volcano\\]", r))); expect_false(any(grepl("!\\[forest\\]", r)))
})

test_that("lme4 is optional and the install hint says how to get it", {
  d <- utils::packageDescription("metamicrobiomeR")
  expect_false(grepl("lme4", d$Imports)); expect_true(grepl("lme4", d$Suggests))
  expect_error(metamicrobiomeR:::.need("notARealPackage123", "a subject random effect"), "Package 'notARealPackage123' is needed for a subject random effect. Install it with install.packages")
})
