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
