# Tests for the importers family; merged from several smaller files (see git history).

# ================ importers ================
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

# ================ importer edge cases ================
cnt3 <- function() matrix(c(10, 5, 3, 2, 8, 6, 4, 9, 7, 1, 12, 5), 3, 4, dimnames = list(c("a", "b", "c"), paste0("S", 1:4)))
smp4 <- function() data.frame(group = factor(c("x", "y", "x", "y")), row.names = paste0("S", 1:4))
info <- list(id = "g", feature_type = "asv")

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

# ================ results from other tools ================
# ---------- effects_from ----------
mq <- data.frame(feature = rep(c("g1", "g2", "g3"), 2), metadata = "group", value = rep(c("b", "c"), each = 3), coef = c(.5, -.2, .1, .8, .0, -.4),
                 stderr = c(.1, .2, .1, .2, .2, .1), pval = c(1e-5, .3, .4, 1e-4, 1, 1e-3), name = rep(c("groupb", "groupc"), each = 3), qval = 0.01, N = 40, N.not.0 = 30, stringsAsFactors = FALSE)
test_that("effects_from reads MaAsLin2, ANCOM-BC2, LinDA and generic tables", {
  e <- effects_from(mq, "maaslin2", "S1", effect_scale = "maaslin2_coef")
  expect_s3_class(e, "mm_effects"); expect_equal(nrow(e), 6); expect_setequal(e$term, c("groupb", "groupc")); expect_true(all(e$engine == "maaslin2"))
  expect_equal(e$estimate[e$feature == "g1" & e$term == "groupb"], 0.5); expect_equal(e$se[e$feature == "g3" & e$term == "groupc"], 0.1)
  expect_equal(e$p_value[e$feature == "g2" & e$term == "groupb"], 0.3); expect_equal(unique(e$n), 40L); expect_equal(unique(e$effect_scale), "maaslin2_coef")
  expect_equal(nrow(effects_from(mq, "maaslin2", "S1", effect_scale = "c", term = "groupb")), 3)
  ab <- data.frame(taxon = c("g1", "g2"), `lfc_(Intercept)` = 1, lfc_groupb = c(.4, -.3), `se_(Intercept)` = .1, se_groupb = c(.1, .2), `p_(Intercept)` = .5, p_groupb = c(.001, .2), q_groupb = .1, check.names = FALSE)
  ea <- effects_from(ab, "ancombc2", "S2", effect_scale = "lfc")
  expect_equal(unique(ea$term), "groupb"); expect_equal(ea$estimate, c(.4, -.3)); expect_equal(ea$engine[1], "ancombc2")
  li <- data.frame(baseMean = c(10, 20), log2FoldChange = c(1, -.5), lfcSE = c(.2, .3), stat = 1, pvalue = c(1e-4, .1), padj = .2, row.names = c("g1", "g2"))
  el <- effects_from(li, "linda", "S3", effect_scale = "log2fc", term = "groupb")
  expect_equal(el$feature, c("g1", "g2")); expect_equal(el$se, c(.2, .3)); expect_error(effects_from(li, "linda", "S3", effect_scale = "log2fc"), "give it with 'term'")
  gen <- data.frame(taxon = c("a", "b"), var = "x", beta = c(.2, .3), sd = c(.1, .1), pp = c(.04, .01), nn = c(10, 12))
  eg <- effects_from(gen, "generic", "S4", effect_scale = "s", columns = list(feature = "taxon", term = "var", estimate = "beta", se = "sd", p_value = "pp", n = "nn"))
  expect_equal(eg$estimate, c(.2, .3)); expect_equal(eg$n, c(10L, 12L))
})

test_that("results from other tools pool, combine, and are kept apart from other engines", {
  mq2 <- mq; mq2$coef <- mq2$coef * 0.9; mq2$stderr <- mq2$stderr * 1.1
  b <- bind_effects(effects_from(mq, "maaslin2", "S1", effect_scale = "coef", term = "groupb"), effects_from(mq2, "maaslin2", "S2", effect_scale = "coef", term = "groupb"))
  p <- pool_effects(b, method = "DL", test = "z", min_prop = 0); expect_equal(nrow(p), 3)
  expect_equal(combine_pvalues(b)$k, c(2, 2, 2))
  other <- new_effect_table("S3", "g1", "groupb", 0.2, 0.1, part = "abundance", engine = "lm_asin", effect_scale = "coef")
  expect_error(bind_effects(b, other), "different engines")
  f <- tempfile(fileext = ".tsv"); utils::write.table(mq, f, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_equal(nrow(effects_from(f, "maaslin2", "S1", effect_scale = "coef")), 6)
})

test_that("effects_from validates its input", {
  expect_error(effects_from(mq, "maaslin2", "S1"), "effect_scale")
  expect_error(effects_from(mq, "maaslin2", effect_scale = "c"), "single study ID")
  expect_error(effects_from(mq[, -4], "maaslin2", "S1", effect_scale = "c"), "lacks the column")
  expect_error(effects_from(mq, "maaslin2", "S1", effect_scale = "c", term = "nope"), "No rows")
  expect_error(effects_from(mq, "generic", "S1", effect_scale = "c"), "columns")
  expect_error(effects_from(data.frame(taxon = "a"), "ancombc2", "S1", effect_scale = "c"), "lfc")
  expect_error(effects_from(data.frame(a = 1), "generic", "S1", effect_scale = "c", columns = list(feature = "x", term = "y", estimate = "z", se = "w", p_value = "v")), "lacks the column")
})
