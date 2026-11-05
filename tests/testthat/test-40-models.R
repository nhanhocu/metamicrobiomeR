# Tests for R/40-model-abundance.R, R/41-model-diversity.R

# ================ fit_study ================
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
  e <- suppressWarnings(fit_study(s, ~ group + age))
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
  e <- suppressWarnings(fit_study(make_bezi_study(), ~ group))
  expect_false("Rare" %in% e$feature)
  e2 <- suppressWarnings(fit_study(make_bezi_study(), ~ group, min_prevalence = 0.01))
  expect_true("Rare" %in% attr(e2, "fit_info")$feature)
})

test_that("a longitudinal model adds a subject random intercept", {
  skip_if_not_installed("gamlss")
  s <- make_bezi_study(subject = "person")
  e <- suppressWarnings(fit_study(s, ~ group))
  expect_type(e$converged, "logical")
  # the feature with a real group effect converges; convergence of the null
  # feature is not guaranteed with this engine
  expect_true(e$converged[e$feature == "F1"])
  expect_setequal(unique(e$term), "groupb")
  expect_gt(e$estimate[e$feature == "F1"], 0.3)
  expect_error(fit_study(make_bezi_study(), ~ group, longitudinal = TRUE), "subject")
})

test_that("fit_study validates its input", {
  s <- make_bezi_study()
  expect_error(fit_study(s, group ~ age), "one-sided")
  expect_error(fit_study(s, ~ nope), "not found in the sample data")
  expect_error(fit_study(s, ~ group, rank = "genus"), "rank")
  expect_error(fit_study(list(), ~ group), "mm_study")
})

test_that("a warning reports fits that did not converge", {
  skip_if_not_installed("gamlss")
  s <- make_bezi_study()
  expect_warning(fit_study(s, ~ group, n_cycles = 2), "did not converge")
})

test_that("features unassigned at their own rank are excluded unless asked for", {
  skip_if_not_installed("gamlss")
  s <- make_bezi_study()
  tx <- data.frame(genus = c("Gen1", NA, "Gen3"), rank = "genus", row.names = rownames(mm_relabund(s)))
  s2 <- as_study(relabund = mm_relabund(s), taxonomy = tx, samples = mm_samples(s),
                 study = list(id = "sim", feature_type = "taxon"))
  e <- suppressWarnings(fit_study(s2, ~ group))
  expect_false("F2" %in% e$feature)
  e2 <- suppressWarnings(fit_study(s2, ~ group, exclude_unassigned = FALSE))
  expect_true("F2" %in% e2$feature)
})

test_that("the lm_asin engine returns effects on the arcsine square root scale", {
  s <- make_bezi_study()
  e <- fit_study(s, ~ group + age, engine = "lm_asin")
  expect_true(all(e$engine == "lm_asin" & e$effect_scale == "asin_sqrt_mean"))
  expect_setequal(unique(e$term), c("groupb", "age"))
  g1 <- e[e$feature == "F1" & e$term == "groupb", ]
  expect_gt(g1$estimate, 0.05)
  expect_lt(g1$p_value, 0.05)
  expect_gt(g1$df, 0)                       # cross-sectional: lm with t-based p-values
  s2 <- make_bezi_study(subject = "person")
  e2 <- fit_study(s2, ~ group, engine = "lm_asin")
  expect_true(all(is.na(e2$df)))            # mixed model: normal approximation, no df
  expect_equal(e2$p_value, 2 * stats::pnorm(-abs(e2$estimate / e2$se)))
  other <- suppressWarnings(fit_study(s, ~ group + age)); other$study <- "sim2"   # a different study, bezi engine
  expect_error(bind_effects(e, other), "effect scales")
})

test_that("fitting on several cores gives identical results", {
  skip_on_os("windows")
  s <- make_bezi_study()
  a <- fit_study(s, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0, cores = 1)
  b <- fit_study(s, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0, cores = 2)
  expect_identical(as.data.frame(a), as.data.frame(b)); expect_identical(attr(a, "fit_info"), attr(b, "fit_info"))
  skip_if_not_installed("gamlss")
  g1 <- suppressWarnings(fit_study(s, ~ group, longitudinal = FALSE, cores = 1)); g2 <- suppressWarnings(fit_study(s, ~ group, longitudinal = FALSE, cores = 2))
  expect_equal(as.data.frame(g1), as.data.frame(g2))
  old <- options(metamicrobiomeR.cores = 2); on.exit(options(old))
  expect_identical(as.data.frame(fit_study(s, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0)), as.data.frame(a))   # the option is used
  expect_error(fit_study(s, ~ group, cores = 0), "'cores' must be a single number of at least 1")
  expect_error(fit_study(s, ~ group, cores = "two"), "'cores' must be")
})

# ================ count engine ================
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

# ================ alpha and beta diversity ================
mk_alpha <- function(n = 40, seed = 3, subject = FALSE) {
  set.seed(seed); grp <- factor(rep(c("a", "b"), length.out = n))
  rich <- ifelse(grp == "b", 60, 25)
  counts <- t(vapply(seq_len(n), function(i) { x <- numeric(80); k <- rich[i]; x[seq_len(k)] <- rpois(k, 40) + 1; x }, numeric(80)))
  dimnames(counts) <- list(paste0("S", seq_len(n)), paste0("F", 1:80)); counts <- t(counts)
  smp <- data.frame(group = grp, row.names = colnames(counts)); if (subject) smp$person <- rep(seq_len(n / 2), each = 2)
  as_study(counts = counts, samples = smp, study = list(id = "al", feature_type = "asv", subject = if (subject) "person"))
}

test_that("alpha_diversity computes the standard indices", {
  cnt <- matrix(c(10, 10, 0, 0,   1, 1, 2, 0,   5, 0, 0, 0), 4, 3, dimnames = list(letters[1:4], c("even", "rare", "mono")))
  s <- as_study(counts = cnt, samples = data.frame(g = 1:3, row.names = colnames(cnt)), study = list(id = "a", feature_type = "asv"))
  a <- alpha_diversity(s, c("shannon", "simpson", "inv_simpson", "pielou", "observed", "chao1"))
  expect_equal(a["even", "shannon"], log(2)); expect_equal(a["even", "simpson"], 0.5); expect_equal(a["even", "inv_simpson"], 2)
  expect_equal(a["even", "pielou"], 1); expect_equal(a["mono", "shannon"], 0); expect_true(is.na(a["mono", "pielou"]))
  expect_equal(a[["observed"]], c(2, 3, 1))
  expect_equal(a["rare", "chao1"], 3 + 2 * 1 / (2 * (1 + 1)))                   # S + f1 (f1 - 1) / (2 (f2 + 1)) with f1 = 2, f2 = 1
  rel <- as_study(relabund = sweep(cnt, 2, colSums(cnt), "/"), samples = mm_samples(s), study = list(id = "r", feature_type = "taxon"))
  expect_equal(alpha_diversity(rel, "shannon")$shannon, a$shannon)
  expect_error(alpha_diversity(rel, "observed"), "need counts")
})

test_that("rarefaction gives the requested depth, is reproducible and drops shallow samples", {
  s <- mk_alpha(10); d <- mm_depth(s); target <- floor(min(d))
  r1 <- alpha_diversity(s, "observed", rarefy = target, seed = 1); r2 <- alpha_diversity(s, "observed", rarefy = target, seed = 1)
  expect_equal(r1, r2); expect_equal(nrow(r1), 10)
  expect_true(all(r1$observed <= alpha_diversity(s, "observed")$observed))
  expect_equal(nrow(alpha_diversity(s, "observed", rarefy = floor(median(d)), seed = 1)), sum(d >= floor(median(d))))
  expect_error(alpha_diversity(s, "observed", rarefy = 1e9), "at least")
})

test_that("fit_alpha returns standardised effects equal to a direct linear model", {
  s <- mk_alpha(); e <- fit_alpha(s, ~ group, metrics = c("shannon", "observed"), longitudinal = FALSE)
  expect_s3_class(e, "mm_effects"); expect_true(all(e$part == "overall" & e$effect_scale == "sd_units" & e$engine == "lm_alpha"))
  al <- alpha_diversity(s, c("shannon", "observed"))
  man <- stats::coef(stats::lm(as.numeric(scale(al$observed)) ~ mm_samples(s)$group))[[2]]
  expect_equal(e$estimate[e$feature == "observed"], man, tolerance = 1e-8)
  expect_gt(e$estimate[e$feature == "observed"], 1)                              # group b is much richer (in SD units)
  raw <- fit_alpha(s, ~ group, metrics = "observed", standardise = FALSE, longitudinal = FALSE)
  expect_equal(raw$effect_scale, "alpha_raw"); expect_gt(raw$estimate, 20)
  expect_equal(raw$estimate / stats::sd(al$observed), e$estimate[e$feature == "observed"], tolerance = 1e-8)
  m <- fit_alpha(mk_alpha(subject = TRUE), ~ group, metrics = "shannon")         # mixed model
  expect_true(all(is.na(m$df)) && nrow(m) == 1)
  expect_error(fit_alpha(s, ~ nope), "not found in the sample data"); expect_error(fit_alpha(s, ~ group, longitudinal = TRUE), "subject")
  # the effects can be pooled across studies like any feature
  e2 <- as.data.frame(e); e2$study <- "al2"; class(e2) <- class(e)
  two <- bind_effects(e, e2)
  expect_equal(nrow(pool_effects(two, method = "DL", test = "z")), 2)
})

mk_beta <- function(n = 40, seed = 8, effect = TRUE, between = FALSE) {
  set.seed(seed)
  person <- rep(seq_len(n / 2), each = 2)
  grp <- if (between) factor(ifelse(person <= n / 4, "a", "b")) else factor(rep(c("a", "b"), length.out = n))
  base <- c(rep(1, 10), rep(0.2, 10)); shift <- c(rep(0.2, 10), rep(1, 10))
  cnt <- vapply(seq_len(n), function(i) stats::rmultinom(1, 2000, if (effect && grp[i] == "b") shift else base), numeric(20))
  dimnames(cnt) <- list(paste0("F", 1:20), paste0("S", seq_len(n)))
  as_study(counts = cnt, samples = data.frame(group = grp, person = person, row.names = colnames(cnt)),
           study = list(id = "be", feature_type = "asv"))
}

test_that("fit_beta finds a composition difference and not a random one", {
  skip_if_not_installed("vegan")
  b <- fit_beta(mk_beta(), ~ group, permutations = 199, seed = 1)
  expect_s3_class(b, "mm_beta"); expect_equal(b$term, "group"); expect_gt(b$R2, 0.3); expect_lt(b$p_value, 0.05)
  nul <- fit_beta(mk_beta(effect = FALSE), ~ group, permutations = 199, seed = 1)
  expect_lt(nul$R2, 0.1); expect_gt(nul$p_value, 0.05)
  a <- fit_beta(mk_beta(), ~ group, distance = "aitchison", permutations = 99, seed = 1)
  expect_gt(a$R2, 0.3); expect_equal(a$distance, "aitchison")
  expect_error(fit_beta(mk_beta(), ~ nope), "not found in the sample data")
})

test_that("with a subject, whole subjects are permuted and R2 is unchanged", {
  skip_if_not_installed("vegan"); skip_if_not_installed("permute")
  s <- mk_beta(between = TRUE); s$study$subject <- "person"
  plain <- fit_beta(s, ~ group, permutations = 99, seed = 2, longitudinal = FALSE)
  blocked <- fit_beta(s, ~ group, permutations = 99, seed = 2)
  expect_equal(blocked$R2, plain$R2, tolerance = 1e-10)
  expect_lt(blocked$p_value, 0.05)
  expect_gte(blocked$p_value, 1 / 100)
})

test_that("a term that varies within subjects is refused instead of returning p = 1", {
  skip_if_not_installed("vegan"); skip_if_not_installed("permute")
  s <- mk_beta(between = FALSE); s$study$subject <- "person"                    # group alternates within each person
  expect_error(fit_beta(s, ~ group, permutations = 49), "within subjects")
  expect_s3_class(fit_beta(s, ~ group, permutations = 49, longitudinal = FALSE), "mm_beta")
  expect_error(fit_beta(mk_beta(), ~ group, longitudinal = TRUE), "needs study\\$subject")
})

test_that("pool_beta combines permutation p-values by Fisher or weighted Stouffer", {
  mkb <- function(study, p, n, r2) structure(data.frame(study = study, term = "group", df = 1, R2 = r2, F = 5, p_value = p, n = n, distance = "bray"), class = c("mm_beta", "data.frame"))
  a <- mkb("s1", 0.01, 30, 0.10); b <- mkb("s2", 0.04, 50, 0.06)
  f <- pool_beta(a, b, method = "fisher")
  expect_equal(f$p_value, stats::pchisq(-2 * (log(0.01) + log(0.04)), 4, lower.tail = FALSE))
  s <- pool_beta(a, b)
  w <- sqrt(c(30, 50)); z <- stats::qnorm(c(0.01, 0.04), lower.tail = FALSE)
  expect_equal(s$p_value, stats::pnorm(sum(w * z) / sqrt(sum(w^2)), lower.tail = FALSE))
  expect_equal(s$R2, stats::weighted.mean(c(0.10, 0.06), c(30, 50))); expect_equal(s$k, 2)
  expect_lt(s$p_value, 0.01)                                         # two supportive studies combine to stronger evidence
  bj <- as.data.frame(b); bj$distance <- "jaccard"; class(bj) <- class(a)
  expect_error(pool_beta(a, bj), "same distance")
  expect_error(pool_beta(a, a), "more than once"); expect_error(pool_beta(data.frame(a = 1)), "mm_beta")
})

# ================ rank selection and failed fits ================
cnt3 <- function() matrix(c(10, 5, 3, 2, 8, 6, 4, 9, 7, 1, 12, 5), 3, 4, dimnames = list(c("a", "b", "c"), paste0("S", 1:4)))
smp4 <- function() data.frame(group = factor(c("x", "y", "x", "y")), row.names = paste0("S", 1:4))
info <- list(id = "g", feature_type = "asv")

multi_rank <- function() {
  tax <- data.frame(kingdom = "Bacteria", phylum = c("Firmicutes", "Bacteroidota", "Firmicutes", "Firmicutes", "Bacteroidota"),
                    genus = c(NA, NA, "Blautia", "Roseburia", "Bacteroides"), rank = c("phylum", "phylum", "genus", "genus", "genus"),
                    row.names = c("Fir", "Bac", "Blau", "Rose", "Bact"), stringsAsFactors = FALSE)
  set.seed(1); n <- 24
  g <- rbeta(n, 4, 6); h <- runif(n, 0.05, 0.3); k <- 1 - g - h
  rel <- rbind(Fir = g * 0 + 0.6, Bac = 0.4, Blau = g * 0.5, Rose = h * 0.5, Bact = k * 0.4)
  colnames(rel) <- paste0("S", 1:n)
  as_study(relabund = rel, taxonomy = tax, samples = data.frame(grp = factor(rep(c("a", "b"), 12)), row.names = colnames(rel)),
           study = list(id = "mr", feature_type = "taxon"))
}

test_that("fit_study, fit_beta and subset functions can select a rank; clr works on proportions", {
  s <- multi_rank()
  e <- fit_study(s, ~ grp, engine = "lm_asin", longitudinal = FALSE, rank = "genus", min_prevalence = 0)
  expect_setequal(unique(e$feature), c("Blau", "Rose", "Bact")); expect_true(all(e$rank == "genus"))
  expect_error(fit_study(s, ~ grp, rank = "family"), "No features with rank")
  plain <- as_study(counts = cnt3(), samples = smp4(), study = info)
  expect_error(fit_study(plain, ~ group, rank = "genus"), "no 'rank' column")
  cl <- fit_study(s, ~ grp, engine = "lm_clr", longitudinal = FALSE, rank = "genus", min_prevalence = 0)         # proportions, no counts
  expect_true(all(cl$effect_scale == "clr_mean")); expect_equal(nrow(cl), 3)
  skip_if_not_installed("vegan")
  b <- fit_beta(s, ~ grp, rank = "genus", permutations = 49, seed = 1)
  expect_equal(b$term, "grp")
  expect_error(fit_beta(s, ~ grp, rank = "family"), "No features with rank")
  expect_error(fit_beta(subset_samples(s, c(1, 2, 3)), ~ grp), "Too few usable samples")
  expect_error(fit_beta(list(), ~ grp), "mm_study"); expect_error(fit_beta(s, "grp"), "one-sided")
  expect_equal(colnames(mm_relabund(subset_samples(s, c(1, 3)))), c("S1", "S3"))                             # numeric selection
  expect_error(subset_samples(s, c(TRUE, FALSE)), "logical vector of length"); expect_error(subset_samples(list(), TRUE), "mm_study")
  copy <- function(x, id) { y <- x; y$study <- id; y }
  three <- bind_effects(list(e, copy(e, "o2"), copy(e, "o3")))
  l <- leave_one_out(three, feature = "Blau", method = "DL", test = "z")
  expect_equal(nrow(l), 3); expect_setequal(l$omitted, c("mr", "o2", "o3")); expect_true(all(l$feature == "Blau"))
  expect_error(leave_one_out(bind_effects(list(e, copy(e, "o2"))), method = "DL", test = "z"), "No feature could be pooled")
})

test_that("failed fits are reported with their error and do not stop the analysis", {
  cn <- cnt3(); one <- data.frame(group = factor(rep("x", 4)), row.names = paste0("S", 1:4))
  s <- as_study(counts = cn, samples = one, study = info)
  e <- expect_warning(fit_study(s, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0), "fits failed")
  expect_equal(nrow(e), 0); expect_s3_class(e, "mm_effects"); expect_false(any(is.na(attr(e, "fit_info")$feature)))
  expect_true(all(!is.na(attr(e, "fit_info")$error)))
  skip_if_not_installed("gamlss")
  e2 <- suppressWarnings(fit_study(s, ~ group, engine = "bezi_gamlss", longitudinal = FALSE, min_prevalence = 0))
  expect_equal(nrow(e2), 0)
  skip_if_not_installed("glmmTMB")
  e3 <- suppressWarnings(fit_study(s, ~ group, engine = "zibb_glmmtmb", longitudinal = FALSE, min_prevalence = 0))
  expect_equal(nrow(e3), 0)
})

# ---- the batch linear fit must give exactly what lm() gives, feature by feature
batch_study <- function(n = 60, nf = 40, na = FALSE) {
  set.seed(3)
  cn <- matrix(rnbinom(nf * n, mu = 10, size = 1), nf, n, dimnames = list(paste0("f", 1:nf), paste0("S", 1:n)))
  cn[2, ] <- c(rep(0L, n - 1), 3L)                                       # an almost empty feature
  smp <- data.frame(group = factor(rep(c("a", "b"), n / 2)), age = rnorm(n), sex = factor(sample(c("F", "M"), n, TRUE)), row.names = colnames(cn))
  if (na) smp$age[c(4, 9)] <- NA
  rel <- sweep(cn, 2, colSums(cn), "/") * 0.9
  rel[1, ] <- 0.05                                                         # exactly constant: a perfect fit, for which lm() warns
  as_study(relabund = rel, samples = smp, study = list(id = "x", feature_type = "taxon"))
}
fit_both <- function(s, formula, engine) {
  old <- options(metamicrobiomeR.batch_lm = FALSE); on.exit(options(old))
  ref <- suppressWarnings(fit_study(s, formula, engine = engine, longitudinal = FALSE, min_prevalence = 0))
  options(metamicrobiomeR.batch_lm = TRUE)
  list(ref = ref, new = suppressWarnings(fit_study(s, formula, engine = engine, longitudinal = FALSE, min_prevalence = 0)))
}

test_that("the batch linear fit equals the per-feature fits for both linear engines", {
  s <- batch_study()
  for (eng in c("lm_asin", "lm_clr")) {
    r <- fit_both(s, ~ group + age + sex, eng)
    expect_equal(as.data.frame(r$new), as.data.frame(r$ref), tolerance = 1e-8)
    expect_equal(attr(r$new, "fit_info"), attr(r$ref, "fit_info"))
  }
  # the helper is really used, and flags the perfectly fitting feature so that lm() handles it as before
  d0 <- mm_samples(s); Y <- asin(sqrt(mm_relabund(s)))
  b <- metamicrobiomeR:::.lm_batch(~ group + age + sex, d0, Y)
  expect_false(is.null(b)); expect_false(b$ok[1]); expect_gte(sum(b$ok), nrow(Y) - 2)
})

test_that("the batch fit steps aside when the design or the data do not allow it", {
  expect_null(metamicrobiomeR:::.lm_batch(~ group, mm_samples(batch_study()), { y <- asin(sqrt(mm_relabund(batch_study()))); y[3, 5] <- NA; y }))   # missing values
  d0 <- mm_samples(batch_study()); Y <- asin(sqrt(mm_relabund(batch_study())))
  expect_null(metamicrobiomeR:::.lm_batch(~ age + I(2 * age), d0, Y))                                 # not full rank
  expect_null(metamicrobiomeR:::.lm_batch(~ group, transform(d0, group = factor("a")), Y))            # one level: lm() reports the failure per feature
  r <- fit_both(batch_study(na = TRUE), ~ group + age, "lm_asin")                                      # missing covariate values
  expect_equal(as.data.frame(r$new), as.data.frame(r$ref), tolerance = 1e-8)
  r2 <- fit_both(batch_study(), ~ group + age + I(2 * age), "lm_asin")                                 # rank-deficient formula
  expect_equal(as.data.frame(r2$new), as.data.frame(r2$ref), tolerance = 1e-8)
})

test_that("when fits fail, the hint says that the rare features fail, with the user's own numbers", {
  h <- metamicrobiomeR:::.convergence_hint(c(a = .06, b = .08, c = .10, d = .12, e = .15, f = .20), c(g = .6, h = .8, i = .9), 0.05)
  expect_match(h, "^ The features that did not converge are rarer \\(present in a median of 11% of the samples, against 80% for the others\\); try a higher 'min_prevalence', for example 0.15\\.$")
  expect_equal(metamicrobiomeR:::.convergence_hint(c(.1, .2, .3), c(.8, .9), 0.05), "")                           # too few failures to say anything
  expect_equal(metamicrobiomeR:::.convergence_hint(rep(.8, 6), rep(.5, 5), 0.05), "")                            # failures are not the rare features
  expect_match(metamicrobiomeR:::.convergence_hint(rep(.7, 6), rep(.95, 3), 0.05), "for example 0.5\\.")          # never suggests more than 0.5
  expect_match(metamicrobiomeR:::.convergence_hint(rep(.05, 6), rep(.9, 3), 0.1), "for example 0.2\\.")           # at least twice the current filter
  expect_match(metamicrobiomeR:::.convergence_hint(rep(.45, 6), rep(.9, 3), 0.05), "for example 0.45\\.")         # otherwise the 75th percentile of the failures, rounded up to 0.05
})

test_that("bezi_gamlss survives a sparse feature whose coefficients are aliased", {
  skip_if_not_installed("gamlss"); skip_on_cran()
  set.seed(3)
  sam <- data.frame(group = rep(c("a", "b"), each = 9), food = rep(c("x", "y", "z"), 6), subj = rep(paste0("s", 1:6), each = 3),
                    row.names = paste0("S", 1:18))
  m <- rbind(common = stats::rbeta(18, 2, 20), rare = c(0.001, rep(0, 8), 0.002, rep(0, 8)))
  m <- rbind(m, other = 1 - colSums(m)); colnames(m) <- rownames(sam)
  st <- as_study(relabund = m, samples = sam, study = list(id = "sparse", feature_type = "otu", subject = "subj"))
  f <- suppressWarnings(suppressMessages(fit_study(st, ~ group + food, longitudinal = TRUE)))
  expect_true(all(c("common", "rare") %in% f$feature))
  expect_false(any(grepl("Intercept", f$term)))
})
