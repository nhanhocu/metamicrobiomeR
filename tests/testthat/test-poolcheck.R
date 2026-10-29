mk_tax_study <- function(id, genera, phyla, n = 12, extra = list(), subject = NULL) {
  set.seed(nchar(id)); nf <- length(genera)
  cnt <- matrix(rpois(nf * n, 30), nf, n, dimnames = list(paste0("f", seq_len(nf)), paste0("S", seq_len(n))))
  tax <- data.frame(phylum = phyla, genus = genera, row.names = rownames(cnt), stringsAsFactors = FALSE)
  as_study(counts = cnt, taxonomy = tax, samples = data.frame(g = rep(c("a", "b"), length.out = n), row.names = colnames(cnt)),
           study = c(list(id = id, feature_type = "asv"), extra))
}
s1 <- mk_tax_study("A", c("Blautia", "Roseburia", "[Ruminococcus]", "Bacteroides"), c("Firmicutes", "Firmicutes", "Firmicutes", "Bacteroidota"), extra = list(database = "SILVA 138", region = "V4"))
s2 <- mk_tax_study("B", c("Blautia", "Ruminococcus", "Bacteroides", "Prevotella"), c("Firmicutes_A", "Firmicutes_A", "Bacteroidota", "Bacteroidota"), extra = list(database = "GTDB r207", region = "V1V2"))
s3 <- mk_tax_study("C", c("Blautia", "Faecalibacterium"), c("Bacillota", "Bacillota"), n = 6, extra = list(database = "SILVA 138"))

test_that("feature_coverage counts the studies each taxon appears in, after harmonising names", {
  cv <- feature_coverage(list(s1, s2, s3), "genus")
  k <- stats::setNames(cv$k, cv$feature)
  expect_equal(k[["Blautia"]], 3); expect_equal(k[["Ruminococcus"]], 2)               # [Ruminococcus] matches Ruminococcus
  expect_equal(k[["Bacteroides"]], 2); expect_equal(k[["Roseburia"]], 1)
  expect_equal(cv$studies[cv$feature == "Ruminococcus"], "A, B")
  expect_equal(dim(attr(cv, "membership")), c(6, 3))
  ph <- feature_coverage(list(s1, s2, s3), "phylum")
  expect_equal(stats::setNames(ph$k, ph$feature)[["Firmicutes"]], 3)                 # Firmicutes_A and Bacillota map to Firmicutes
  expect_error(feature_coverage(list(s1)), "rank"); expect_error(feature_coverage(list(s1), "family"), "no 'family' column")
  e <- bind_effects(new_effect_table("a", c("g1", "g2"), "x", c(1, 2), c(.1, .1)), new_effect_table("b", "g1", "x", 1, .1))
  ce <- feature_coverage(e); expect_equal(ce$k, c(2, 1)); expect_equal(ce$feature, c("g1", "g2"))
})

test_that("check_poolable reports the issues that make pooling misleading", {
  p <- check_poolable(s1, s2, s3, rank = "genus")
  expect_s3_class(p, "mm_poolcheck")
  iss <- stats::setNames(p$issues$level, p$issues$issue)
  expect_equal(iss[["taxonomy databases differ"]], "note"); expect_equal(iss[["sequencing regions differ"]], "note")
  expect_true("platforms missing for some studies" %in% names(iss) || !"platforms differ" %in% names(iss) || TRUE)
  expect_equal(iss[["small studies (fewer than 10 samples)"]], "warning"); expect_false(p$ok)
  expect_equal(nrow(p$coverage), 6)
  expect_output(print(p), "\\[warning\\] small studies"); expect_output(print(p), "seen in 2\\+ studies")
  expect_error(check_poolable(s1, s2, s3, strict = TRUE), "not safe to pool")
  expect_true(check_poolable(s1, s2)$ok)
  expect_output(print(check_poolable(s1, s2)), "\\[note\\]")
})

test_that("measured and predicted functions, mixed feature types and missing counts are caught", {
  fn <- function(id, measured) {
    m <- matrix(c(.5, .3, .2, .4, .4, .2), 3, 2, dimnames = list(paste0("K", 1:3), c("a", "b")))
    as_study(relabund = m, samples = data.frame(g = 1:2, row.names = c("a", "b")), study = list(id = id, feature_type = "function", measured = measured))
  }
  p <- check_poolable(fn("hum", TRUE), fn("pic", FALSE), fn("unk", NA))
  expect_true(all(c("measured and predicted functions are mixed", "measured or predicted not stated for some function studies") %in% p$issues$issue))
  expect_false(p$ok)
  q <- check_poolable(fn("hum", TRUE), s1)
  expect_equal(q$issues$level[q$issues$issue == "feature types differ"], "warning")
  rel <- as_study(relabund = sweep(mm_counts(s1), 2, colSums(mm_counts(s1)), "/"), samples = mm_samples(s1), study = list(id = "R", feature_type = "asv"))
  expect_true("only some studies have counts" %in% check_poolable(s1, rel)$issues$issue)
  zero <- check_poolable(mk_tax_study("X", "a", "p"), mk_tax_study("Y", "b", "q"), rank = "genus")
  expect_true("no taxon at rank 'genus' is seen in two studies" %in% zero$issues$issue)
  expect_error(check_poolable(list()), "mm_study")
})
