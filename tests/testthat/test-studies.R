# Tests for the studies family; merged from several smaller files (see git history).

# ================ study object ================
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

# ================ malformed input and hand-edited studies ================
cnt3 <- function() matrix(c(10, 5, 3, 2, 8, 6, 4, 9, 7, 1, 12, 5), 3, 4, dimnames = list(c("a", "b", "c"), paste0("S", 1:4)))
smp4 <- function() data.frame(group = factor(c("x", "y", "x", "y")), row.names = paste0("S", 1:4))
info <- list(id = "g", feature_type = "asv")

test_that("as_study accepts counts and proportions together and aligns them", {
  cn <- cnt3(); rel <- sweep(cn, 2, colSums(cn), "/")
  s <- as_study(counts = as.data.frame(cn), relabund = rel[c("c", "a", "b"), c("S4", "S3", "S2", "S1")], samples = smp4(), study = info)
  expect_equal(names(s$assays), c("counts", "relabund"))
  expect_identical(dimnames(s$assays$relabund), dimnames(s$assays$counts))
  expect_equal(unname(mm_relabund(s)["a", ]), unname(rel["a", ]))                       # reordered to match the counts
  expect_error(as_study(counts = cn, relabund = rel[1:2, ], samples = smp4(), study = info), "IDs do not match between 'counts' and 'relabund'")
})

test_that("as_study rejects malformed input", {
  cn <- cnt3()
  expect_error(as_study(counts = "x", samples = smp4(), study = info), "numeric matrix")
  expect_error(as_study(counts = cn, samples = list(), study = info), "data.frame")
  expect_error(as_study(counts = cn, samples = smp4(), study = 1), "named list")
  dup <- cn; rownames(dup)[2] <- "a"
  expect_error(as_study(counts = dup, samples = smp4(), study = info), "unique")
  expect_error(as_study(counts = cn, taxonomy = "x", samples = smp4(), study = info), "data.frame")
  expect_error(as_study(counts = cn, taxonomy = data.frame(phylum = "p", row.names = "zzz"), samples = smp4(), study = info), "do not match")
  na <- cn; na[1, 1] <- NA
  expect_error(as_study(counts = na, samples = smp4(), study = info), "missing value")
})

test_that("check_study catches a study that was modified by hand", {
  s <- as_study(counts = cnt3(), samples = smp4(), study = info)
  bad <- s; bad$assays <- list(); expect_error(check_study(bad), "at least one assay")
  bad <- s; bad$assays$relabund <- bad$assays$counts[, 4:1]; expect_error(check_study(bad), "identical dimnames")
  bad <- s; rownames(bad$features)[1] <- "zzz"; expect_error(check_study(bad), "features' rows")
  bad <- s; rownames(bad$samples)[1] <- "zzz"; expect_error(check_study(bad), "samples' rows")
  expect_error(check_study(list()), "mm_study"); expect_error(study_table(), "mm_study")
})


# ================ merge_studies ================
# ---------- merge_studies ----------
mk_tax <- function(id, genera, phyla, n = 12, subject = FALSE, counts = TRUE) {
  set.seed(nchar(id) + 7); nf <- length(genera)
  cnt <- matrix(rpois(nf * n, 30), nf, n, dimnames = list(paste0("f", seq_len(nf)), paste0("S", seq_len(n))))
  smp <- data.frame(group = factor(rep(c("a", "b"), length.out = n)), row.names = colnames(cnt))
  if (subject) smp$person <- rep(seq_len(n / 2), each = 2)
  tax <- data.frame(phylum = phyla, genus = genera, row.names = rownames(cnt), stringsAsFactors = FALSE)
  as_study(counts = if (counts) cnt else NULL, relabund = if (!counts) sweep(cnt, 2, colSums(cnt), "/") else NULL, taxonomy = tax, samples = smp,
           study = list(id = id, feature_type = "asv", subject = if (subject) "person"))
}
sa <- mk_tax("A", c("Blautia", "Roseburia", "[Ruminococcus]", "Bacteroides"), c("Firmicutes", "Firmicutes", "Firmicutes", "Bacteroidota"), subject = TRUE)
sb <- mk_tax("B", c("Blautia", "Ruminococcus", "Bacteroides", "Prevotella"), c("Firmicutes_A", "Firmicutes_A", "Bacteroidota", "Bacteroidota"), n = 10, subject = TRUE)

test_that("merge_studies matches taxa across databases and builds one study", {
  m <- merge_studies(list(sa, sb), rank = "genus")
  expect_s3_class(m, "mm_study")
  expect_setequal(rownames(mm_counts(m)), c("Blautia", "Bacteroides", "Ruminococcus"))                  # intersection by harmonised name
  expect_equal(ncol(mm_counts(m)), 22); expect_true(all(c("A:S1", "B:S10") %in% colnames(mm_counts(m))))
  expect_equal(levels(mm_samples(m)$study), c("A", "B")); expect_equal(as.character(mm_samples(m)["B:S3", "study"]), "B")
  # [Ruminococcus] and Ruminococcus are a single taxon within and across studies
  expect_equal(unname(mm_counts(m)["Ruminococcus", "A:S2"]), unname(mm_counts(sa)["f3", "S2"]))
  expect_equal(unname(mm_counts(m)["Blautia", "B:S4"]), unname(mm_counts(sb)["f1", "S4"]))
  expect_equal(mm_study_info(m)$id, "A+B"); expect_equal(mm_study_info(m)$merged_from, c("A", "B"))
  expect_equal(mm_study_info(m)$subject, "subject_merged"); expect_equal(mm_samples(m)["A:S1", "subject_merged"], "A:1")
  expect_equal(length(unique(mm_samples(m)$subject_merged)), 6 + 5)                                      # subjects are not shared across studies
  expect_true(all(mm_features(m)$rank == "genus")); expect_equal(mm_features(m)["Blautia", "genus"], "Blautia")
})

test_that("intersection and union differ in how missing features are handled", {
  u <- merge_studies(list(sa, sb), rank = "genus", features = "union")
  expect_setequal(rownames(mm_counts(u)), c("Blautia", "Roseburia", "Ruminococcus", "Bacteroides", "Prevotella"))
  expect_true(all(mm_counts(u)["Prevotella", grep("^A:", colnames(mm_counts(u)))] == 0))                 # absent from A's database: filled with zero
  expect_true(all(mm_counts(u)["Roseburia", grep("^B:", colnames(mm_counts(u)))] == 0))
  expect_equal(unname(colSums(mm_counts(u))[grep("^A:", colnames(mm_counts(u)))]), unname(colSums(mm_counts(sa))))   # totals preserved
})

test_that("merge_studies handles proportions, missing subjects and invalid input", {
  rel <- mk_tax("R", c("Blautia", "Bacteroides", "Roseburia"), c("F", "B", "F"), counts = FALSE)
  expect_warning(m <- merge_studies(list(sa, rel), rank = "genus"), "merging relative abundances")
  expect_null(mm_counts(m)); expect_true(!is.null(mm_relabund(m)))
  nosub <- mk_tax("N", c("Blautia", "Bacteroides"), c("F", "B"), subject = FALSE)
  expect_warning(merge_studies(list(sa, nosub), rank = "genus"), "subject column is missing")
  expect_error(merge_studies(list(sa), rank = "genus"), "at least two")
  expect_error(merge_studies(list(sa, sa), rank = "genus"), "unique")
  other <- mk_tax("O", c("Zzz"), "P"); expect_error(merge_studies(list(sa, other), rank = "genus"), "No feature is shared")
  clash <- sb; clash$samples$study <- "x"; expect_error(merge_studies(list(sa, clash), rank = "genus"), "already a sample column")
})

test_that("one-stage and two-stage analyses agree on simulated studies with a common effect", {
  sim <- simulate_studies(n_studies = 3, tau = 0, n_samples = 80, n_features = 12, n_diff = 3, effect = 1.5, depth_sd = 0.3, seed = 11)
  eff <- lapply(sim$studies, fit_study, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0)
  two <- pool_effects(bind_effects(eff), term = "groupb", method = "DL", test = "z", min_prop = 0)
  one <- fit_study(merge_studies(sim$studies), ~ group + study, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0)
  one <- one[one$term == "groupb", ]
  j <- merge(two, one, by = "feature")
  expect_gt(stats::cor(j$estimate.x, j$estimate.y), 0.9)
  expect_equal(sign(j$estimate.x[abs(j$estimate.x) > 0.05]), sign(j$estimate.y[abs(j$estimate.x) > 0.05]))
  dif <- sim$truth$feature[sim$truth$diff]
  expect_true(all(one$p_value[one$feature %in% dif] < 0.01)); expect_equal(unique(one$study), "sim1+sim2+sim3")
})

test_that("merge_studies keeps the order of factor levels, so the reference level and the sign of effects do not change", {
  mk <- function(id, seed) {
    set.seed(seed); cn <- matrix(rpois(6 * 40, 30), 6, 40, dimnames = list(paste0("g", 1:6), paste0("S", 1:40)))
    grp <- factor(rep(c("control", "case"), each = 20), levels = c("control", "case")); cn[1, grp == "case"] <- cn[1, grp == "case"] * 3L
    as_study(counts = cn, taxonomy = data.frame(genus = paste0("g", 1:6), row.names = rownames(cn)), samples = data.frame(group = grp, row.names = colnames(cn)),
             study = list(id = id, feature_type = "taxon"))
  }
  s <- list(mk("A", 1), mk("B", 2))
  m <- suppressMessages(merge_studies(s, rank = "genus"))
  expect_equal(levels(mm_samples(m)$group), c("control", "case"))                                  # not re-sorted alphabetically ("case" first)
  per <- bind_effects(lapply(s, fit_study, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0))
  one <- suppressWarnings(fit_study(m, ~ group + study, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0))
  expect_true("groupcase" %in% per$term); expect_true("groupcase" %in% one$term); expect_false("groupcontrol" %in% one$term)
  a <- per$estimate[per$term == "groupcase" & per$feature == "g1" & per$study == "A"]; b <- one$estimate[one$term == "groupcase" & one$feature == "g1"]
  expect_equal(sign(a), sign(b)); expect_gt(b, 0)
  mixed <- list(s[[1]], { x <- s[[2]]; x$samples$group <- factor(c("case", "control")[as.integer(x$samples$group)], levels = c("case", "control")); x })   # different order in one study: first study wins
  expect_equal(levels(mm_samples(suppressMessages(merge_studies(mixed, rank = "genus")))$group), c("control", "case"))
  chr <- lapply(s, function(x) { x$samples$group <- as.character(x$samples$group); x })                                                              # plain text columns still become factors
  expect_s3_class(mm_samples(suppressMessages(merge_studies(chr, rank = "genus")))$group, "factor")
})
