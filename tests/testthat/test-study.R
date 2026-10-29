make_toy <- function(n_feat = 6, n_samp = 8, seed = 1) {
  set.seed(seed)
  counts <- matrix(rpois(n_feat * n_samp, 20), n_feat, n_samp,
                   dimnames = list(paste0("F", seq_len(n_feat)), paste0("S", seq_len(n_samp))))
  tax <- data.frame(Phylum = rep(c("A", "B"), each = n_feat / 2),
                    Genus = c("g1", "g2", "g3", "g4", NA, ""),
                    row.names = rownames(counts))
  samp <- data.frame(group = rep(c("x", "y"), n_samp / 2), subj = rep(1:4, each = 2),
                     row.names = colnames(counts))
  list(counts = counts, tax = tax, samp = samp)
}
info <- list(id = "toy", feature_type = "asv", region = "V4")

test_that("as_study builds a valid object from counts", {
  d <- make_toy()
  s <- as_study(counts = d$counts, taxonomy = d$tax, samples = d$samp, study = c(info, subject = "subj"))
  expect_s3_class(s, "mm_study")
  expect_equal(dim(mm_counts(s)), c(6, 8))
  expect_equal(unname(colSums(mm_relabund(s))), rep(1, 8))
  expect_equal(unname(mm_depth(s)), unname(colSums(d$counts)))
  expect_output(print(s), "<mm_study> toy")
})

test_that("taxonomy rank columns are normalised and rank is never inferred", {
  d <- make_toy()
  s <- as_study(counts = d$counts, taxonomy = d$tax, samples = d$samp, study = info)
  f <- mm_features(s)
  expect_true(all(c("phylum", "genus") %in% colnames(f)))
  expect_false("rank" %in% colnames(f))
  expect_true(is.na(f$genus[6]))  # empty string becomes NA
})

test_that("a user-supplied rank groups features for the closure check", {
  # phylum-level and genus-level features, each level closed separately
  rel <- rbind(P1 = c(0.6, 0.5), P2 = c(0.4, 0.5), G1 = c(0.3, 0.2), G2 = c(0.7, 0.8))
  colnames(rel) <- c("s1", "s2")
  smp <- data.frame(g = 1:2, row.names = colnames(rel))
  tx <- data.frame(rank = c("phylum", "phylum", "genus", "genus"), row.names = rownames(rel))
  inf <- list(id = "ml", feature_type = "taxon")
  expect_s3_class(as_study(relabund = rel, taxonomy = tx, samples = smp, study = inf), "mm_study")
  expect_error(as_study(relabund = rel, samples = smp, study = inf), "sum to at most 1")
})

test_that("samples and taxonomy are aligned to the assay order", {
  d <- make_toy()
  s <- as_study(counts = d$counts, taxonomy = d$tax[6:1, ], samples = d$samp[8:1, ], study = info)
  expect_identical(rownames(mm_features(s)), rownames(mm_counts(s)))
  expect_identical(rownames(mm_samples(s)), colnames(mm_counts(s)))
})

test_that("percent input is converted to proportions", {
  rel <- matrix(c(60, 40, 25, 75), 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  s <- as_study(relabund = rel, samples = data.frame(g = 1:2, row.names = c("s1", "s2")),
                study = list(id = "pct", feature_type = "taxon"), relabund_scale = "percent")
  expect_equal(unname(colSums(mm_relabund(s))), c(1, 1))
  expect_true(all(is.na(mm_depth(s))))
})

test_that("proportions that do not close within a rank are rejected, partial sums are fine", {
  bad <- matrix(c(60, 40, 25, 75), 2, dimnames = list(c("a", "b"), c("s1", "s2")))
  smp <- data.frame(g = 1:2, row.names = c("s1", "s2"))
  expect_error(as_study(relabund = bad, samples = smp, study = list(id = "x", feature_type = "taxon")),
               "percent")
  ok <- bad / 100 * 0.5
  expect_s3_class(as_study(relabund = ok, samples = smp, study = list(id = "x", feature_type = "taxon")), "mm_study")
})

test_that("invalid input is rejected with clear errors", {
  d <- make_toy()
  expect_error(as_study(samples = d$samp, study = info), "counts")
  expect_error(as_study(counts = d$counts, samples = d$samp, study = list(id = "x")), "feature_type")
  expect_error(as_study(counts = d$counts, samples = d$samp, study = list(feature_type = "asv")), "id")
  expect_error(as_study(counts = d$counts, samples = d$samp, study = list(id = "x", feature_type = "banana")), "feature_type")
  neg <- d$counts; neg[1, 1] <- -1
  expect_error(as_study(counts = neg, samples = d$samp, study = info), "negative")
  frac <- d$counts + 0.5
  expect_error(as_study(counts = frac, samples = d$samp, study = info), "whole numbers")
  expect_error(as_study(counts = d$counts, samples = d$samp[-1, ], study = info), "samples")
  nn <- unname(d$counts)
  expect_error(as_study(counts = nn, samples = d$samp, study = info), "names")
  expect_error(as_study(counts = d$counts, samples = d$samp, study = c(info, subject = "nope")), "subject")
})

test_that("zero-depth samples give a warning, not an error", {
  d <- make_toy(); d$counts[, 3] <- 0L
  expect_warning(s <- as_study(counts = d$counts, samples = d$samp, study = info), "zero total")
  expect_true(all(is.na(mm_relabund(s)[, 3])))
})

test_that("NULL study metadata falls back to the defaults instead of breaking the check", {
  d <- make_toy()
  s <- as_study(counts = d$counts, samples = d$samp,
                study = list(id = "toy", feature_type = "asv", subject = NULL, region = NULL))
  expect_true(is.na(mm_study_info(s)$subject))
  expect_true(is.na(mm_study_info(s)$region))
})
