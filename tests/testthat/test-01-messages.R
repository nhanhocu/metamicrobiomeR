# Tests for R/01-messages.R
reset_once <- function() assign("told", character(), envir = metamicrobiomeR:::.mm_state)
sim1 <- function(...) simulate_study(n_samples = 40, n_features = 8, seed = 1, ...)$study
fit_q <- function(s, ...) fit_study(s, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0, ...)

# ------------------------------------------------------------ helpers
test_that("the helpers format messages consistently", {
  expect_equal(metamicrobiomeR:::.show(letters[1:3]), "a, b, c")
  expect_equal(metamicrobiomeR:::.show(letters[1:7]), "a, b, c, ... (4 more)")
  expect_equal(metamicrobiomeR:::.plural(1, "sample"), "1 sample"); expect_equal(metamicrobiomeR:::.plural(2, "study", "studies"), "2 studies")
  expect_equal(metamicrobiomeR:::.pct(0.5), "50%")
  expect_error(metamicrobiomeR:::.check_class(1, "mm_study", "study", "as_study()"), "'study' must be an mm_study object \\(see as_study\\(\\)\\); got a numeric")
  expect_error(metamicrobiomeR:::.check_class(NULL, "data.frame", "x"), "must be a data.frame object; got a NULL")
  expect_equal(metamicrobiomeR:::.id_mismatch("Sample IDs", "A", c("a", "b"), "B", c("b", "c")), "Sample IDs do not match between A and B: 1 only in A (a); 1 only in B (c).")
  expect_null(metamicrobiomeR:::.id_mismatch("x", "A", "a", "B", "a"))
  expect_error(metamicrobiomeR:::.need("notARealPackage123", "testing"), "Install it with install.packages\\(\"notARealPackage123\"\\)")
  expect_error(metamicrobiomeR:::.missing_vars("zz", c("a", "b")), "Variable not found in the sample data: zz. Available: a, b.")
  expect_error(metamicrobiomeR:::.check_list_of(list(1), "mm_study", "as_study()"), "Supply mm_study objects \\(see as_study\\(\\)\\); got a numeric")
})

test_that("messages can be silenced with an option, warnings cannot, and repeated messages appear once", {
  old <- options(metamicrobiomeR.quiet = NULL); on.exit(options(old))
  na <- sim1(); na$samples$group[1:6] <- NA
  expect_message(fit_q(na), "6 of 40 samples have missing values")
  options(metamicrobiomeR.quiet = TRUE)
  expect_silent(fit_q(na))
  expect_warning(fit_study(sim1(), ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 5), "No feature passed the filters")   # still warns
  options(metamicrobiomeR.quiet = NULL); reset_once()
  sub <- simulate_study(n_samples = 40, n_features = 8, n_subjects = 20, between = TRUE, seed = 1)$study
  expect_message(suppressMessages(fit_study(sub, ~ group, engine = "lm_asin")), NA)
  reset_once()
  m1 <- testthat::capture_messages(fit_study(sub, ~ group, engine = "lm_asin")); m2 <- testthat::capture_messages(fit_study(sub, ~ group, engine = "lm_asin"))
  expect_true(any(grepl("Using 'person' as a subject random intercept", m1))); expect_false(any(grepl("Using 'person'", m2)))      # once per session
})

# ------------------------------------------------------------ study object
test_that("study-object errors say which IDs and which values are the problem", {
  cn <- matrix(1:6, 2, 3, dimnames = list(c("a", "b"), c("S1", "S2", "S3")))
  smp <- data.frame(g = 1:3, row.names = c("S1", "S2", "S9"))
  expect_error(as_study(counts = cn, samples = smp, study = list(id = "x", feature_type = "asv")), "Sample IDs do not match between the abundance table and 'samples': 1 only in the abundance table \\(S3\\); 1 only in 'samples' \\(S9\\)")
  dup <- cn; rownames(dup)[2] <- "a"
  expect_error(as_study(counts = dup, samples = smp[1:3, , drop = FALSE], study = list(id = "x", feature_type = "asv")), "Feature IDs must be unique; 'a' is repeated")
  rel <- sweep(cn, 2, colSums(cn), "/") * 50
  expect_error(as_study(relabund = rel, samples = data.frame(g = 1:3, row.names = colnames(cn)), study = list(id = "x", feature_type = "taxon")), "the largest sum is 50. Percent data\\?")
  expect_error(as_study(counts = cn + 0.5, samples = data.frame(g = 1:3, row.names = colnames(cn)), study = list(id = "x", feature_type = "asv")), "whole numbers \\(6 values are not\\)")
  expect_error(as_study(counts = cn, samples = data.frame(g = 1:3, row.names = colnames(cn)), study = list(id = "x", feature_type = "nonsense")), "got 'nonsense'")
  expect_error(as_study(counts = cn, samples = data.frame(g = 1:3, row.names = colnames(cn)), study = list(id = "x", feature_type = "asv", subject = "zz")), "Available: g")
  expect_error(as_study(counts = cn, samples = data.frame(g = 1:3, row.names = colnames(cn)), study = list(feature_type = "asv")), "study\\$id is required; for example")
  z <- cn; z[, 2] <- 0L
  expect_warning(as_study(counts = z, samples = data.frame(g = 1:3, row.names = colnames(cn)), study = list(id = "x", feature_type = "asv")), "1 sample has zero total abundance \\(S2\\); kept, but left out of models")
  expect_error(as_study(counts = "x", samples = smp, study = list(id = "x", feature_type = "asv")), "got a character")
})

test_that("a 'Domain' column is the kingdom, and a taxonomy without ranks is flagged", {
  cn <- matrix(1:6, 2, 3, dimnames = list(c("a", "b"), c("S1", "S2", "S3"))); smp <- data.frame(g = 1:3, row.names = colnames(cn))
  s <- as_study(counts = cn, taxonomy = data.frame(Domain = "Bacteria", Phylum = c("P", "Q"), row.names = c("a", "b")), samples = smp, study = list(id = "x", feature_type = "asv"))
  expect_true(all(c("kingdom", "phylum") %in% colnames(mm_features(s))))
  expect_warning(as_study(counts = cn, taxonomy = data.frame(foo = c("p", "q"), row.names = c("a", "b")), samples = smp, study = list(id = "x", feature_type = "asv")), "has no rank columns")
  expect_silent(as_study(counts = cn, taxonomy = data.frame(sequence = c("A", "C"), row.names = c("a", "b")), samples = smp, study = list(id = "x", feature_type = "asv")))   # not meant as taxonomy
})

# ------------------------------------------------------------ fitting
test_that("fit_study reports what it filtered, what it left out, and why nothing was fitted", {
  s <- sim1()
  expect_warning(e <- fit_study(s, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 2, min_mean = 1), "No feature passed the filters \\(present in at least 200% of the samples")
  expect_equal(nrow(e), 0)
  expect_message(fit_study(s, ~ group, engine = "lm_asin", longitudinal = FALSE, min_mean = 0.2), "Fitting \\d of 8 features \\(\\d below the prevalence and abundance filters\\)")
  expect_silent(fit_q(s))                                                                             # nothing notable: nothing printed
  expect_error(fit_study(s, ~ grp), "Variable not found in the sample data: grp. Available: group.")
  expect_error(fit_study(mm_counts(s), ~ group), "'study' must be an mm_study object \\(see as_study\\(\\)\\); got a matrix")
  expect_error(fit_study(s, ~ group, rank = "genus"), "no 'rank' column, so 'rank' cannot select features; agglomerate\\(\\) creates one")
})

test_that("singular mixed-model fits are reported once, not once per feature", {
  sub <- simulate_study(n_samples = 40, n_features = 8, n_subjects = 20, subject_sd = 0.001, between = TRUE, seed = 1)$study
  reset_once()
  msgs <- testthat::capture_messages(e <- fit_study(sub, ~ group, engine = "lm_asin", min_prevalence = 0))
  expect_lte(length(msgs), 2)
  expect_true(any(grepl("mixed-model fits were singular", msgs)))
  expect_true("singular" %in% colnames(attr(e, "fit_info"))); expect_true("first_warning" %in% colnames(attr(e, "fit_info")))
})

test_that("failed and non-converged fits give a short warning with the next step", {
  cn <- matrix(c(10, 5, 3, 2, 8, 6, 4, 9, 7, 1, 12, 5), 3, 4, dimnames = list(c("a", "b", "c"), paste0("S", 1:4)))
  one <- as_study(counts = cn, samples = data.frame(group = factor(rep("x", 4)), row.names = colnames(cn)), study = list(id = "o", feature_type = "asv"))
  expect_warning(fit_study(one, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0), "3 of 3 fits failed \\(first error: .*\\); see attr\\(, \"fit_info\"\\)")
  skip_if_not_installed("gamlss")
  expect_warning(fit_study(sim1(), ~ group, engine = "bezi_gamlss", longitudinal = FALSE, n_cycles = 2, min_prevalence = 0), "did not converge within 2 iterations; try a larger 'n_cycles'")
})

# ------------------------------------------------------------ pooling
mkq <- function(study, est, se, feature = "g", engine = "t", scale = "s") new_effect_table(study, feature, "x", est, se, rank = "genus", part = "abundance", n = 20L, engine = engine, effect_scale = scale)
test_that("pooling refuses estimates that cannot be compared, even without bind_effects", {
  mixed <- rbind(as.data.frame(mkq("a", 0.5, 0.1)), as.data.frame(mkq("b", 0.4, 0.1, engine = "other"))); class(mixed) <- class(mkq("a", 1, 1))
  expect_error(pool_effects(mixed), "Studies were fitted with different engines for abundance / x \\(other vs t\\). Fit every study with the same engine")
  sc <- rbind(as.data.frame(mkq("a", 0.5, 0.1)), as.data.frame(mkq("b", 0.4, 0.1, scale = "z"))); class(sc) <- class(mixed)
  expect_error(pool_effects(sc), "different effect scales for abundance / x \\(s vs z\\)")
  expect_error(meta_regress(sc, data.frame(study = c("a", "b"), r = c("u", "v")), ~ r), "different effect scales")
})

test_that("pooling says which features were not pooled, and why nothing could be", {
  b <- bind_effects(mkq(paste0("s", 1:4), c(.5, .6, .4, .5), rep(.1, 4), "g1"), mkq("s1", .2, .1, "g2"))
  expect_message(p <- pool_effects(b, method = "DL", test = "z"), "1 feature of 2 not pooled \\(estimated in fewer than 2 studies or in fewer than 50% of the studies reporting the term\\)")
  expect_equal(p$feature, "g1")
  only <- bind_effects(mkq("s1", .2, .1, "g2"), mkq("s2", .3, .1, "g3"))
  expect_error(pool_effects(only), "No feature was estimated in enough studies to pool: it needs at least 2 studies and 50% of the studies reporting the term, but the most any feature has is 1. Lower 'min_studies' or 'min_prop'")
  two <- bind_effects(mkq(c("a", "b"), c(.5, .6), c(.1, .1)))
  reset_once()
  expect_message(pool_effects(two), "1 feature was pooled from fewer than 5 studies \\(2 studies\\).*with 2 studies it has 1 degree of freedom.*compare with test = \"z\" or run pool_sensitivity\\(\\)")
  expect_silent(suppressMessages(pool_effects(two)))
  expect_silent(pool_effects(two, test = "z"))
  reset_once()
  four <- bind_effects(mkq(paste0("s", 1:4), c(.5, .6, .4, .5), rep(.1, 4), "g1"), mkq(paste0("s", 1:3), c(.5, .6, .4), rep(.1, 3), "g2"))
  m <- testthat::capture_messages(pool_effects(four))
  expect_match(paste(m, collapse = ""), "2 features were pooled from fewer than 5 studies \\(3 to 4 studies\\)")
  expect_false(grepl("1 degree of freedom", paste(m, collapse = "")))
  reset_once()
  five <- bind_effects(mkq(paste0("s", 1:5), c(.5, .6, .4, .5, .55), rep(.1, 5), "g1"))
  expect_false(any(grepl("fewer than 5 studies", testthat::capture_messages(pool_effects(five)))))
  w <- suppressWarnings(bind_effects(mkq(c("a", "b", "c"), c(.5, NA, .4), c(.1, .1, 0))))
  expect_warning(suppressMessages(pool_effects(w, min_studies = 1, min_prop = 0)), "2 estimates not used: 1 with a missing estimate or standard error, 1 with a standard error that is not positive")
  expect_warning(new_effect_table("s", "f", "x", 1, 0), "1 standard error is exactly zero; such estimates cannot be weighted")
  expect_error(bind_effects(mkq("a", .5, .1), mkq("a", .4, .1)), "Duplicate rows in the effect table: study 'a', feature 'g', part 'abundance', term 'x' appears more than once")
  expect_error(pool_effects(two, term = "nope"), "No effects left after selecting term 'nope'")
})

test_that("effect and pooled tables print a short header instead of a raw data frame", {
  b <- bind_effects(mkq(paste0("s", 1:4), c(.5, .6, .4, .5), rep(.1, 4), "g1"))
  out <- testthat::capture_output(print(b))
  expect_match(out, "<mm_effects> 4 studies \\(s1, s2, s3, s4\\), 1 feature, 4 rows"); expect_match(out, "engine: t \\| scale: s \\| parts: abundance \\| terms: x")
  expect_match(out, "converged: 4 of 4")
  expect_s3_class(b[1:2, ], "mm_effects")                                                          # subsetting keeps the class
  expect_match(testthat::capture_output(print(b[, c("feature", "estimate")])), "feature estimate")          # dropped columns: no crash, plain table
  p0 <- pool_effects(b, method = "DL", test = "z")
  expect_match(testthat::capture_output(print(p0[, c("feature", "estimate")])), "feature estimate")
  expect_match(testthat::capture_output(print(p0[0, ])), "feature")
  big <- bind_effects(mkq(paste0("s", 1:8), rep(.5, 8), rep(.1, 8), "g1")); expect_match(testthat::capture_output(print(big)), "# ... 2 more rows")
  p <- pool_effects(b, method = "DL", test = "z")
  expect_match(testthat::capture_output(print(p)), "<mm_pooled> 1 feature pooled \\(DL, z test\\); 4 studies each")
  expect_match(testthat::capture_output(print(p)), "adjusted p <= 0.05: 1 \\| median I2: 0%")
  dom <- bind_effects(mkq(c("a", "b", "c"), c(.4, .0, .5), c(.005, .4, .35)))
  expect_match(testthat::capture_output(print(suppressWarnings(pool_effects(dom, method = "DL", test = "z")))), "1 dominated by one study")
})

# ------------------------------------------------------------ other functions
test_that("agglomeration, merging, rarefaction and beta diversity say what they left out", {
  cn <- matrix(rpois(30, 30), 5, 6, dimnames = list(letters[1:5], paste0("S", 1:6)))
  tax <- data.frame(genus = c("A", "B", NA, "C", NA), row.names = letters[1:5])
  s <- as_study(counts = cn, taxonomy = tax, samples = data.frame(g = 1:6, row.names = colnames(cn)), study = list(id = "t", feature_type = "asv"))
  expect_message(agglomerate(s, "genus"), "2 of 5 features \\(\\d+% of the abundance\\) have no assignment at 'genus'; kept as 'unassigned_genus' and left out of models")
  expect_message(agglomerate(s, "genus", unassigned = "drop"), "Dropped 2 of 5 features")
  expect_error(agglomerate(s, "family"), "has no 'family' column \\(it has: genus\\)")
  expect_error(agglomerate(agglomerate(suppressMessages(s), "genus"), "genus"), "already agglomerated \\(its features are at rank genus\\)")
  a <- simulate_study(n_samples = 20, n_features = 10, seed = 1)$study; b <- simulate_study(n_samples = 20, n_features = 10, seed = 2, id = "sim2")$study
  expect_message(merge_studies(list(a, b)), "Merged 2 studies: 10 features \\(of 10, those present in all studies\\) in 40 samples")
  expect_error(merge_studies(list(a)), "needs at least two studies; got 1")
  half <- subset_features(b, as.integer(sub("F", "", rownames(mm_features(b)))) <= 3)
  expect_warning(suppressMessages(merge_studies(list(a, half))), "Only 3 of 10 features \\(30%\\) are present in all studies; consider features = \"union\"")
  d <- mm_depth(a)
  expect_message(alpha_diversity(a, "observed", rarefy = floor(median(d)), seed = 1), "Rarefying to .* reads: \\d+ of 20 samples with fewer reads are dropped")
  expect_error(alpha_diversity(a, "observed", rarefy = 1e9), "No sample has at least 1e\\+09 reads; the deepest has")
  skip_if_not_installed("vegan")
  na <- a; na$samples$group[1:3] <- NA
  expect_message(fit_beta(na, ~ group, permutations = 19, seed = 1), "3 of 20 samples are left out \\(missing model variables or no abundance\\)")
})

test_that("importers report what they could not match", {
  seqtab <- matrix(c(5, 6, 7, 8), 2, 2, dimnames = list(c("S1", "S2"), c("AAA", "CCC")))
  smp <- data.frame(g = 1:2, row.names = c("S1", "S2"))
  taxa <- matrix(c("B", "B", "x", "y"), 2, 2, dimnames = list(c("GGG", "TTT"), c("Kingdom", "Genus")))
  expect_warning(read_dada2(seqtab, taxa, smp, list(id = "d")), "row names of 'taxa' do not match the sequences of 'seqtab'; assuming the rows are in the same order")
  expect_error(read_dada2(seqtab, taxa[1, , drop = FALSE], smp, list(id = "d")), "'taxa' must have one row per sequence of 'seqtab' \\(1 rows for 2 sequences\\)")
  d <- tempfile(); dir.create(d)
  writeLines(c("#OTU ID\tS1\tS2", "f1\t5\t6", "f2\t1\t2", "f3\t3\t3"), file.path(d, "ft.tsv")); writeLines(c("Feature ID\tTaxon", "f1\td__Bacteria;p__Firmicutes"), file.path(d, "tax.tsv"))
  expect_message(read_qiime2(file.path(d, "ft.tsv"), file.path(d, "tax.tsv"), smp, list(id = "q")), "2 of 3 features have no row in the taxonomy table and are left unassigned")
  pf <- data.frame(clade_name = c("k__B|p__F|c__C|o__O|f__F|g__G|s__G_a", "k__B|p__F"), S1 = c(60, 100), S2 = c(30, 100), check.names = FALSE)
  expect_message(read_metaphlan(pf, smp, list(id = "m")), "Kept 1 of 2 rows: the clades at level 'species'")
  expect_error(read_metaphlan(pf, smp, list(id = "m"), level = "genus"), "No clades at level 'genus' in the profile \\(deepest levels found:")
  reset_once()
  ft <- data.frame(id = c("K1", "K2"), S1 = c(1, 3), S2 = c(2, 2))
  expect_message(read_function_table(ft, smp, list(id = "f"), measured = TRUE), "normalised to proportions")
  expect_silent(read_function_table(ft, smp, list(id = "f"), measured = TRUE))                          # once per session
  expect_error(read_function_table(ft, smp, list(id = "f"), measured = "yes"), "measured \\(TRUE\\) or predicted \\(FALSE\\)")
  expect_error(read_qiime2("no-such-file", NULL, smp, list(id = "q")), "File not found: no-such-file")
})

test_that("filtering, hierarchical testing and meta-regression report what they dropped or skipped", {
  eff <- bind_effects(mkq(paste0("s", 1:3), c(.5, .6, .4), rep(.1, 3), "ok"), mkq(paste0("s", 1:3), c(.5, .6, .4), c(.1, .1, .0001), "bad"))
  expect_message(f <- filter_effects(eff), "Dropped 1 of 6 estimates before pooling: 1 with an implausibly small standard error \\(see attr\\(, \"dropped\"\\)\\)")
  expect_silent(filter_effects(mkq(c("a", "b"), c(.5, .6), c(.1, .1))))
  pl <- structure(data.frame(feature = c("P", "C"), rank = c("phylum", "class"), part = "abundance", term = "x", k = 3, estimate = .5, se = .1, ci_lo = 0, ci_hi = 1, statistic = 1,
                             p_value = c(0.5, 1e-6), tau2 = 0, I2 = 0, Q = 0, Q_p = 1, pi_lo = NA_real_, pi_hi = NA_real_, max_weight = .4, top_study = "a", p_adjust = NA_real_), class = c("mm_pooled", "data.frame"))
  expect_message(adjust_hierarchical(pl, c(P = NA, C = "P")), "1 of 2 features were tested; 1 were not because their parent was not rejected. 0 rejected at q = 0.05")
  sd_ <- data.frame(study = paste0("s", 1:3), r = c("u", "u", "u"))
  expect_error(meta_regress(eff, sd_, ~ r), "needs more studies than model coefficients and a moderator that varies between them")
  mixed <- bind_effects(mkq(paste0("s", 1:4), c(.5, .6, .4, .5), rep(.1, 4), "g1"), mkq(c("s1", "s2"), c(.5, .6), c(.1, .1), "g2"))
  expect_message(meta_regress(mixed, data.frame(study = paste0("s", 1:4), r = c("u", "u", "v", "v")), ~ r, test = "z"), "1 of 2 features skipped")
})

test_that("notes are not repeated per setting and never mention a 0% rule", {
  b <- bind_effects(mkq(paste0("s", 1:4), c(.5, .6, .4, .5), rep(.1, 4), "g1"), mkq("s1", .2, .1, "g2"))
  m <- testthat::capture_messages(pool_sensitivity(b, min_prop = 0.5))
  expect_equal(sum(grepl("not pooled", m)), 1)                                           # four settings, one note
  m0 <- testthat::capture_messages(pool_effects(b, min_prop = 0, min_studies = 2))
  expect_match(m0, "1 feature of 2 not pooled \\(estimated in fewer than 2 studies\\)\\.")
  expect_false(any(grepl("0%", m0)))
  res <- suppressMessages(meta_analyse(simulate_studies(n_studies = 3, n_samples = 14, n_features = 6, n_diff = 1, effect = 0.2, tau = 0.5, seed = 2)$studies, ~ group, term = "groupb",
                                       engine = "lm_asin", fit_args = list(longitudinal = FALSE, min_prevalence = 0)))
  expect_match(testthat::capture_output(print(res)), "pass the correction \\(see pool_sensitivity")
})
