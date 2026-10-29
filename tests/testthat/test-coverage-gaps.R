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

test_that("read_metaphlan reads a file with comment lines; read_dada2 works without taxonomy", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("#mpa_vJan21_CHOCOPhlAnSGB", "clade_name\tS1\tS2", "k__Bacteria|p__Firmicutes|c__C|o__O|f__F|g__G|s__G_a\t60\t30",
               "k__Bacteria|p__Bacteroidota|c__C|o__O|f__F|g__H|s__H_b\t40\t70"), f)
  s <- read_metaphlan(f, data.frame(g = 1:2, row.names = c("S1", "S2")), list(id = "mf"))
  expect_equal(unname(mm_relabund(s)[, "S1"]), c(0.6, 0.4)); expect_equal(mm_study_info(s)$feature_type, "taxon")
  writeLines(c("a\tb", "1\t2"), f)
  expect_error(read_metaphlan(f, data.frame(g = 1, row.names = "b"), list(id = "x")), "clade_name")
  expect_error(read_metaphlan("no-such-file", data.frame(), list(id = "x")), "File not found")
  seqtab <- matrix(c(5, 6, 7, 8), 2, 2, dimnames = list(c("S1", "S2"), c("AAA", "CCC")))
  s2 <- read_dada2(seqtab, taxa = NULL, samples = data.frame(g = 1:2, row.names = c("S1", "S2")), study = list(id = "d"))
  expect_equal(rownames(mm_counts(s2)), c("ASV1", "ASV2")); expect_equal(mm_features(s2)$sequence, c("AAA", "CCC"))
  s3 <- read_dada2(seqtab, taxa = cbind(Kingdom = c("B", "B"), Genus = c("x", "y")), samples = data.frame(g = 1:2, row.names = c("S1", "S2")), study = list(id = "d"))
  expect_equal(mm_features(s3)$genus, c("x", "y"))                                              # taxa given by position
})

test_that("a phyloseq object with samples in rows converts", {
  skip_if_not_installed("phyloseq")
  m <- cnt3(); ps <- phyloseq::phyloseq(phyloseq::otu_table(t(m), taxa_are_rows = FALSE))
  s <- study_from_phyloseq(ps, list(id = "p"))
  expect_equal(unname(mm_counts(s)), unname(m)); expect_equal(nrow(mm_samples(s)), 4)
})

test_that("effect-table validation reports what is wrong", {
  e <- new_effect_table("s", c("f1", "f2"), "x", c(0.1, 0.2), c(0.1, 0.1))
  expect_error(check_effects(data.frame(a = 1)), "mm_effects")
  bad <- as.data.frame(e)[, -1]; class(bad) <- class(e); expect_error(check_effects(bad), "missing from the effect table")
  bad <- as.data.frame(e); bad$study[1] <- NA; class(bad) <- class(e); expect_error(check_effects(bad), "must not be missing")
  bad <- as.data.frame(e); bad$converged[1] <- NA; class(bad) <- class(e); expect_error(check_effects(bad), "converged")
  na <- bind_effects(new_effect_table("s1", "f", "x", NA, 0.1), new_effect_table("s2", "f", "x", NA, 0.1))
  expect_error(suppressWarnings(pool_effects(na)), "No usable estimates")
})
