smp3 <- function(ids) data.frame(group = rep(c("a", "b"), length.out = length(ids)), row.names = ids)

test_that("read_dada2 builds a counts study with taxonomy and sequences", {
  set.seed(1)
  seqs <- c("ACGTAC", "ACGTTT", "GGGTAC", "TTTTAC")
  seqtab <- matrix(rpois(12, 40), 3, 4, dimnames = list(c("S1", "S2", "S3"), seqs))
  taxa <- cbind(Kingdom = "Bacteria", Phylum = c("Firmicutes", "Firmicutes", "Bacteroidota", NA),
                Genus = c("Blautia", NA, "Bacteroides", NA)); rownames(taxa) <- rev(seqs)[c(4, 3, 2, 1)]
  rownames(taxa) <- seqs
  s <- read_dada2(seqtab, taxa, smp3(rownames(seqtab)), list(id = "d", region = "V4"))
  expect_equal(dim(mm_counts(s)), c(4, 3)); expect_equal(rownames(mm_counts(s)), paste0("ASV", 1:4))
  expect_equal(unname(mm_counts(s)["ASV2", ]), unname(seqtab[, "ACGTTT"]))
  f <- mm_features(s)
  expect_equal(f$sequence, seqs); expect_equal(f$genus, c("Blautia", NA, "Bacteroides", NA)); expect_equal(mm_study_info(s)$feature_type, "asv")
  expect_error(read_dada2(seqtab, taxa[1:2, ], smp3(rownames(seqtab)), list(id = "d")), "one row per sequence")
  expect_error(read_dada2(unname(seqtab), taxa, smp3(rownames(seqtab)), list(id = "d")), "sample names")
})

test_that("read_qiime2 reads the exported feature table, taxonomy and metadata", {
  d <- tempfile(); dir.create(d)
  writeLines(c("# Constructed from biom file", "#OTU ID\tS1\tS2\tS3", "f1\t10\t20\t30", "f2\t5\t0\t7", "f3\t1\t1\t0"), file.path(d, "ft.tsv"))
  writeLines(c("Feature ID\tTaxon\tConfidence", "f1\td__Bacteria;p__Firmicutes;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia\t0.99",
               "f2\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__;f__;g__\t0.8", "f3\tUnassigned\t0.5"), file.path(d, "tax.tsv"))
  writeLines(c("sample-id\tgroup", "#q2:types\tcategorical", "S1\ta", "S2\tb", "S3\ta"), file.path(d, "meta.tsv"))
  s <- read_qiime2(file.path(d, "ft.tsv"), file.path(d, "tax.tsv"), file.path(d, "meta.tsv"), list(id = "q"))
  expect_equal(dim(mm_counts(s)), c(3, 3)); expect_equal(unname(mm_counts(s)["f2", ]), c(5, 0, 7))
  f <- mm_features(s)
  expect_equal(f["f1", "genus"], "Blautia"); expect_equal(f["f2", "phylum"], "Bacteroidota"); expect_true(is.na(f["f2", "genus"]))
  expect_true(all(is.na(f["f3", c("kingdom", "phylum")])))
  expect_equal(mm_samples(s)$group, c("a", "b", "a")); expect_equal(rownames(mm_samples(s)), c("S1", "S2", "S3"))
  expect_error(read_qiime2(file.path(d, "ft.tsv"), data.frame(a = 1), smp3(c("S1", "S2", "S3")), list(id = "q")), "Feature ID")
})

test_that("read_metaphlan keeps one level and converts percent to proportions", {
  pf <- data.frame(clade_name = c("k__Bacteria", "k__Bacteria|p__Firmicutes", "k__Bacteria|p__Firmicutes|c__Clostridia|o__E|f__L|g__Blautia|s__Blautia_obeum",
                                  "k__Bacteria|p__Firmicutes|c__Clostridia|o__E|f__L|g__Blautia|s__Blautia_obeum|t__SGB1",
                                  "k__Bacteria|p__Bacteroidota|c__B|o__B|f__B|g__Bacteroides|s__Bacteroides_fragilis"),
                   S1 = c(100, 60, 40, 40, 55), S2 = c(100, 70, 20, 20, 75), check.names = FALSE)
  s <- read_metaphlan(pf, smp3(c("S1", "S2")), list(id = "m"), level = "species")
  expect_equal(nrow(mm_relabund(s)), 2)                                   # strain-level row excluded
  expect_equal(unname(mm_relabund(s)[1, ]), c(0.40, 0.20))
  expect_true(all(mm_features(s)$rank == "species")); expect_equal(mm_features(s)$genus[2], "Bacteroides")
  expect_error(read_metaphlan(pf, smp3(c("S1", "S2")), list(id = "m"), level = "family"), "No clades")
})

test_that("phyloseq and SummarizedExperiment objects convert and round-trip", {
  skip_if_not_installed("phyloseq"); skip_if_not_installed("SummarizedExperiment")
  data("GlobalPatterns", package = "phyloseq", envir = environment())
  gp <- phyloseq::prune_taxa(phyloseq::taxa_names(GlobalPatterns)[1:60], phyloseq::prune_samples(phyloseq::sample_names(GlobalPatterns)[1:6], GlobalPatterns))
  s <- study_from_phyloseq(gp, list(id = "gp", region = "V4"))
  expect_equal(dim(mm_counts(s)), c(60, 6))
  expect_true(all(c("kingdom", "phylum", "genus") %in% colnames(mm_features(s))))
  expect_equal(unname(mm_counts(s)[1, ]), unname(as(phyloseq::otu_table(gp), "matrix")[1, ]))
  expect_true("SampleType" %in% colnames(mm_samples(s)))
  back <- as_phyloseq(s)
  expect_equal(unname(as(phyloseq::otu_table(back), "matrix")), unname(mm_counts(s)))
  se <- as_summarized_experiment(s)
  s2 <- study_from_se(se, "counts", list(id = "se"))
  expect_equal(unname(mm_counts(s2)), unname(mm_counts(s))); expect_equal(mm_features(s2)$phylum, mm_features(s)$phylum)
  rel <- phyloseq::transform_sample_counts(gp, function(x) x / sum(x))
  expect_null(mm_counts(study_from_phyloseq(rel, list(id = "r"))))
})

test_that("read_biom_study reads a BIOM file", {
  skip_if_not_installed("biomformat")
  m <- matrix(c(5, 0, 3, 2, 8, 1), 3, 2, dimnames = list(c("o1", "o2", "o3"), c("S1", "S2")))
  obs <- data.frame(taxonomy1 = "k__Bacteria", taxonomy2 = c("p__Firmicutes", "p__Bacteroidetes", "p__Firmicutes"), row.names = rownames(m))
  b <- biomformat::make_biom(m, observation_metadata = obs)
  f <- tempfile(fileext = ".biom"); biomformat::write_biom(b, f)
  s <- read_biom_study(f, smp3(c("S1", "S2")), list(id = "b"))
  expect_equal(unname(mm_counts(s)["o3", ]), c(3, 1)); expect_equal(mm_features(s)["o2", "phylum"], "Bacteroidetes")
  expect_equal(mm_study_info(s)$feature_type, "otu")
})

test_that("read_function_table needs an explicit measured flag and normalises samples", {
  d <- data.frame(`function` = c("K00001", "K00002", "K00003"), description = c("a", "b", "c"), S1 = c(10, 30, 60), S2 = c(1, 1, 2), check.names = FALSE)
  sm <- data.frame(nsti = c(0.05, 0.4), row.names = c("S1", "S2"))
  expect_error(read_function_table(d, sm, list(id = "f")), "measured")
  s <- read_function_table(d, sm, list(id = "f"), measured = FALSE, source = "PICRUSt2 2.5.2")
  expect_equal(unname(colSums(mm_relabund(s))), c(1, 1)); expect_equal(unname(mm_relabund(s)["K00003", ]), c(0.6, 0.5))
  i <- mm_study_info(s); expect_equal(i$feature_type, "function"); expect_identical(i$measured, FALSE); expect_equal(i$pipeline, "PICRUSt2 2.5.2")
  expect_equal(mm_features(s)$description, c("a", "b", "c"))
  low <- subset_samples(s, nsti < 0.1)
  expect_equal(colnames(mm_relabund(low)), "S1"); expect_equal(nrow(mm_samples(low)), 1)
  expect_error(subset_samples(s, nsti < 0), "No samples")
})

test_that("subset_samples and subset_features work on counts and proportions", {
  set.seed(2); cnt <- matrix(rpois(20, 20), 5, 4, dimnames = list(letters[1:5], paste0("S", 1:4)))
  s <- as_study(counts = cnt, taxonomy = data.frame(phylum = c("A", "A", "B", "B", "C"), row.names = letters[1:5]),
                samples = data.frame(g = c("x", "y", "x", "y"), row.names = colnames(cnt)), study = list(id = "s", feature_type = "asv"))
  a <- subset_samples(s, g == "x"); expect_equal(colnames(mm_counts(a)), c("S1", "S3"))
  b <- subset_features(s, phylum %in% c("A", "B")); expect_equal(rownames(mm_counts(b)), letters[1:4])
  rel <- as_study(relabund = mm_relabund(s), samples = mm_samples(s), study = list(id = "r", feature_type = "taxon"))
  expect_equal(colnames(mm_relabund(subset_samples(rel, g == "y"))), c("S2", "S4"))
  expect_equal(colnames(mm_relabund(subset_samples(s, c(TRUE, FALSE, FALSE, FALSE)))), "S1")
  expect_error(subset_samples(s, g == "z"), "No samples"); expect_error(subset_features(list(), TRUE), "mm_study")
})
