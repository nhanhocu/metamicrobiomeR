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

# ---------- overview plots ----------
mkp <- function(study, est, se, feature, p = NA_real_) new_effect_table(study = study, feature = feature, term = "x", estimate = est, se = se, rank = "genus", part = "abundance",
                                                                        p_value = p, n = 20L, engine = "t", effect_scale = "logit_mean")
fe <- paste0("g", 1:6)
set.seed(3)
effp <- bind_effects(mkp("A", c(0.9, 0.5, 0.1, -0.2, 0.0, 0.3), rep(0.1, 6), fe, c(1e-9, 1e-4, 0.4, 0.04, 1, 0.01)), mkp("B", c(0.8, 0.6, 0.0, -0.1, 0.1), rep(0.12, 5), fe[1:5], c(1e-8, 1e-3, 1, 0.3, 0.5)),
                     mkp("C", c(1.0, 0.4, 0.2, -0.3, 0.2, 0.2), rep(0.15, 6), fe, c(1e-6, 0.01, 0.2, 0.05, 0.3, 0.2)))
geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
layer <- function(p, geom, nth = 1) ggplot2::ggplot_build(p)$data[[which(geoms(p) == geom)[nth]]]

test_that("plot_heatmap shows every feature in every study plus the pooled column", {
  pl <- pool_effects(effp, term = "x", min_prop = 0)
  h <- plot_heatmap(effp, pl)
  expect_s3_class(h, "ggplot")
  tiles <- layer(h, "GeomTile"); expect_equal(nrow(tiles), 6 * 4)                          # 6 features x (3 studies + pooled)
  expect_equal(sum(is.na(tiles$fill) | tiles$fill == "grey88"), 1)                         # g6 was not estimated in study B
  lab <- layer(h, "GeomText")$label
  expect_true(all(c("***", "*", "") %in% lab)); expect_gt(sum(lab == "***"), 2)
  ord <- levels(h$data$y); expect_equal(rev(ord)[1], pl$feature[which.min(pl$p_value)])    # ordered by pooled p-value
  expect_equal(length(levels(plot_heatmap(effp, pl, top = 3)$data$y)), 3)
  expect_equal(levels(plot_heatmap(effp, pl, features = c("g2", "g4"))$data$y), c("g4", "g2"))
  expect_equal(nrow(layer(plot_heatmap(effp, pl, stars = FALSE), "GeomTile")), 24); expect_false("GeomText" %in% geoms(plot_heatmap(effp, pl, stars = FALSE)))
  expect_equal(levels(plot_heatmap(effp, pl, features = "g1", feature_labels = c(g1 = "Blautia"))$data$y), "Blautia")
  dup <- plot_heatmap(effp, pl, features = c("g1", "g2", "g3"), feature_labels = c(g1 = "same", g2 = "same", g3 = "other"))   # equal labels must not break the plot
  expect_equal(levels(dup$data$y), c("other", "same.1", "same"))
})

test_that("plot_heatmap and plot_volcano validate input and use the pooled table", {
  expect_error(plot_heatmap(data.frame(a = 1)), "mm_effects")
  expect_error(plot_heatmap(effp, features = "zzz"), "not found in the effect table")
  expect_error(plot_heatmap(effp, pooled = data.frame(a = 1)), "mm_pooled")
  expect_s3_class(plot_heatmap(effp, features = "g1"), "ggplot")                            # pooled computed when not supplied
  pl <- pool_effects(effp, term = "x", min_prop = 0)
  v <- plot_volcano(pl, alpha = 0.05, label_top = 2)
  expect_s3_class(v, "ggplot")
  pts <- layer(v, "GeomPoint"); expect_equal(nrow(pts), nrow(pl))
  expect_equal(sort(pts$y), sort(-log10(pl$p_value)), tolerance = 1e-8)
  expect_equal(sum(layer(v, "GeomText")$label != ""), 2)
  expect_error(plot_volcano(data.frame(a = 1)), "mm_pooled")
  ext <- bind_effects(effp, mkp("A", -40, 0.1, "zz"), mkp("B", -40, 0.1, "zz"), mkp("C", -40, 0.1, "zz"))        # a degenerate fit with an extreme estimate
  pe <- pool_effects(ext, term = "x", min_prop = 0, method = "DL", test = "z")
  ve <- plot_volcano(pe); pts <- layer(ve, "GeomPoint")
  expect_equal(sum(pts$shape == 17), 1); expect_lt(max(abs(pts$x)), 10)                                           # clamped to the plotted range, marked
  expect_match(ve$labels$caption, "beyond the plotted range")
  all_ <- layer(plot_volcano(pe, x_range = Inf), "GeomPoint"); expect_equal(min(all_$x), pe$estimate[pe$feature == "zz"], tolerance = 1e-8)
  expect_null(plot_volcano(pool_effects(effp, term = "x", min_prop = 0))$labels$caption)
  prev <- as.data.frame(mkp("A", c(0.2, 0.3, 0.1), rep(0.1, 3), fe[1:3])); prev$part <- "prevalence"; class(prev) <- class(effp)
  prev2 <- prev; prev2$study <- "B"
  two <- pool_effects(bind_effects(effp, prev, prev2), min_prop = 0, term = "x")
  expect_equal(length(unique(ggplot2::ggplot_build(plot_volcano(two))$layout$layout$PANEL)), 2)
})
