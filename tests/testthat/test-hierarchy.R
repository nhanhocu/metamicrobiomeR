pooled_tab <- function(feature, p, rank, part = "abundance", term = "x") {
  x <- data.frame(feature = feature, rank = rank, part = part, term = term, k = 3, estimate = 0.5, se = 0.1, ci_lo = 0, ci_hi = 1,
                  statistic = 1, p_value = p, tau2 = 0, I2 = 0, Q = 0, Q_p = 1, pi_lo = NA_real_, pi_hi = NA_real_,
                  max_weight = 0.4, top_study = "s1", p_adjust = NA_real_, stringsAsFactors = FALSE)
  class(x) <- c("mm_pooled", "data.frame"); x
}
feat <- c("P1", "P2", "C1a", "C1b", "C2a", "C2b", "G1", "G2", "G3")
par  <- c(P1 = NA, P2 = NA, C1a = "P1", C1b = "P1", C2a = "P2", C2b = "P2", G1 = "C1a", G2 = "C1a", G3 = "C2a")
pv   <- c(1e-4, 0.6, 0.001, 0.2, 0.0005, 0.3, 1e-5, 0.5, 1e-6)

test_that("hierarchical FDR tests only the children of rejected features, with the scaled level", {
  h <- adjust_hierarchical(pooled_tab(feat, pv, c("phylum", "phylum", rep("class", 4), rep("genus", 3))), par, q = 0.05)
  g <- stats::setNames(h$p_hier, h$feature)
  expect_equal(unname(h$level), c(1, 1, 2, 2, 2, 2, 3, 3, 3))
  # level 1: BH over the two phyla; one of two rejected, so level 2 is scaled by 1/2
  expect_equal(g[["P1"]], 2e-4); expect_equal(g[["P2"]], 0.6)
  expect_equal(g[["C1a"]], 0.002 / 0.5); expect_equal(g[["C1b"]], 0.2 / 0.5)
  expect_true(is.na(g[["C2a"]]) && is.na(g[["C2b"]]))                   # parent P2 not rejected: not tested
  # level 3 is scaled by (1/2) * (1/4): one of 2 phyla and one of 4 classes rejected
  expect_equal(g[["G1"]], 2e-5 / 0.125); expect_equal(g[["G2"]], min(0.5 / 0.125, 1))
  expect_true(is.na(g[["G3"]]))                                          # despite p = 1e-6: its ancestors were not rejected
  expect_equal(h$feature[h$hier_rejected], c("P1", "C1a", "G1"))
  expect_equal(attr(h, "settings")$hierarchical_q, 0.05)
})

test_that("nothing below an unrejected level is tested, and orphans become roots", {
  h <- adjust_hierarchical(pooled_tab(feat[1:4], c(0.5, 0.6, 1e-9, 1e-9), c("phylum", "phylum", "class", "class")), par, q = 0.05)
  expect_equal(sum(h$hier_rejected), 0); expect_true(all(is.na(h$p_hier[3:4])))
  # a class whose phylum is not in the table is treated as a root
  o <- adjust_hierarchical(pooled_tab(c("C1a", "C1b"), c(0.001, 0.2), c("class", "class")), par, q = 0.05)
  expect_equal(o$level, c(1, 1)); expect_equal(o$feature[o$hier_rejected], "C1a")
})

test_that("families are separate for each part and term", {
  a <- pooled_tab(c("P1", "P2"), c(1e-4, 0.6), "phylum", part = "abundance")
  b <- pooled_tab(c("P1", "P2"), c(0.7, 0.8), "phylum", part = "prevalence")
  h <- adjust_hierarchical(rbind(a, b), par)
  expect_equal(h$feature[h$hier_rejected & h$part == "abundance"], "P1"); expect_false(any(h$hier_rejected[h$part == "prevalence"]))
})

test_that("input is validated", {
  expect_error(adjust_hierarchical(data.frame(a = 1), par), "mm_pooled")
  expect_error(adjust_hierarchical(pooled_tab("P1", 0.1, "phylum"), unname(par)), "named")
  expect_error(adjust_hierarchical(pooled_tab("P1", 0.1, "phylum"), par, q = 2), "'q'")
})

test_that("feature_parents finds the nearest higher feature with the same lineage", {
  tax <- data.frame(kingdom = "Bacteria", phylum = c("Firmicutes", "Bacteroidota", "Firmicutes", "Firmicutes", "Bacteroidota", "Firmicutes"),
                    class = c(NA, NA, "Clostridia", "Clostridia", "Bacteroidia", "Bacilli"),
                    genus = c(NA, NA, NA, "Blautia", NA, "Lactobacillus"),
                    rank = c("phylum", "phylum", "class", "genus", "class", "genus"),
                    row.names = c("Fir", "Bac", "Clo", "Blau", "Bcl", "Lacto"), stringsAsFactors = FALSE)
  rel <- matrix(c(0.6, 0.4, 0.5, 0.3, 0.4, 0.2), 6, 2, dimnames = list(rownames(tax), c("S1", "S2")))
  s <- as_study(relabund = rel * 0.5, taxonomy = tax, samples = data.frame(g = 1:2, row.names = c("S1", "S2")), study = list(id = "h", feature_type = "taxon"))
  p <- feature_parents(s)
  expect_true(is.na(p[["Fir"]]) && is.na(p[["Bac"]]))
  expect_equal(p[["Clo"]], "Fir"); expect_equal(p[["Bcl"]], "Bac")
  expect_equal(p[["Blau"]], "Clo")                       # genus: class is the nearest present higher rank
  expect_equal(p[["Lacto"]], "Fir")                      # no Bacilli class feature: falls back to the phylum
})
