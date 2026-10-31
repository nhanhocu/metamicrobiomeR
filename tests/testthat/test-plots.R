# Tests for the plots family; merged from several smaller files (see git history).

# ================ forest plot ================
mk <- function(study, feature, est, se, term = "x", part = "abundance", scale = "logit_mean") {
  new_effect_table(study = study, feature = feature, term = term, estimate = est, se = se, rank = "genus",
                   part = part, n = 30L, engine = "test", effect_scale = scale)
}
eff <- bind_effects(
  mk("Bangladesh", c("g1", "g2"), c(0.50, -0.10), c(0.10, 0.20)),
  mk("Haiti",       c("g1", "g2"), c(0.40,  0.20), c(0.20, 0.30)),
  mk("USA CA/FL",   c("g1"),       c(0.60),        c(0.15)),
  mk("USA NC",      c("g1", "g2"), c(0.30,  0.00), c(0.30, 0.40)))
titles <- function(p) as.character(ggplot2::ggplot_build(p)$layout$layout$panel)
geoms <- function(p) vapply(p$layers, function(l) class(l$geom)[1], character(1))
layer <- function(p, geom, nth = 1) ggplot2::ggplot_build(p)$data[[which(geoms(p) == geom)[nth]]]

test_that("plot_forest returns a ggplot with one panel per feature", {
  p1 <- plot_forest(eff, feature = "g1")
  expect_s3_class(p1, "ggplot")
  expect_length(titles(p1), 1)
  p2 <- plot_forest(eff, feature = c("g1", "g2"))
  expect_length(titles(p2), 2)
  expect_true(grepl("^g1", titles(p2)[1]) && grepl("^g2", titles(p2)[2]))
})

test_that("the study layer holds the estimates and intervals, the pooled layer the diamond", {
  p <- plot_forest(eff, feature = "g1", level = 0.95)
  seg <- layer(p, "GeomSegment", 1)                            # per-study intervals
  expect_equal(nrow(seg), 4)
  lo <- c(0.50, 0.40, 0.60, 0.30) - stats::qnorm(0.975) * c(0.10, 0.20, 0.15, 0.30)
  hi <- c(0.50, 0.40, 0.60, 0.30) + stats::qnorm(0.975) * c(0.10, 0.20, 0.15, 0.30)
  expect_equal(sort(seg$x), sort(lo), tolerance = 1e-8)
  expect_equal(sort(seg$xend), sort(hi), tolerance = 1e-8)
  pooled <- pool_effects(eff, term = "x")
  poly <- layer(p, "GeomPolygon")                              # pooled diamond
  expect_equal(range(poly$x), c(pooled$ci_lo[pooled$feature == "g1"], pooled$ci_hi[pooled$feature == "g1"]), tolerance = 1e-8)
  expect_equal(unique(poly$x[poly$y > 0]), pooled$estimate[pooled$feature == "g1"], tolerance = 1e-8)
})

test_that("a different confidence level changes the study intervals", {
  w <- function(l) { s <- layer(plot_forest(eff, feature = "g1", level = l), "GeomSegment", 1)[1, ]; s$xend - s$x }
  expect_gt(w(0.99), w(0.90))
})

test_that("a supplied pooled table is used and the panel title reports the pooled result", {
  pl <- pool_effects(eff, method = "DL", test = "z")
  ttl <- titles(plot_forest(eff, pooled = pl, feature = "g1"))
  expect_match(ttl, "DL, z test", fixed = TRUE)
  expect_match(ttl, "k = 4", fixed = TRUE)
  expect_match(ttl, "I2 = ", fixed = TRUE)
  est <- trimws(formatC(pl$estimate[pl$feature == "g1"], digits = 2, format = "g"))
  expect_match(ttl, paste0("pooled ", est, " ["), fixed = TRUE)                              # no stray padding between the label and the number
  default <- titles(plot_forest(eff, feature = "g1"))
  expect_match(default, "REML, adhoc test", fixed = TRUE)
})

test_that("the prediction interval is drawn only with three or more studies, and can be switched off", {
  n_seg <- function(p) sum(geoms(p) == "GeomSegment")                              # intervals + prediction interval
  expect_equal(n_seg(plot_forest(eff, feature = "g1")), 2)                         # g1 has 4 studies
  expect_equal(n_seg(plot_forest(eff, feature = "g1", show_pi = FALSE)), 1)
  expect_equal(n_seg(plot_forest(eff, feature = "g2")), 2)                         # g2 has 3 studies
  e2 <- bind_effects(mk("s1", "f", 0.3, 0.1), mk("s2", "f", 0.2, 0.1))
  expect_equal(n_seg(plot_forest(e2, feature = "f")), 1)                           # 2 studies: no prediction interval
})

test_that("weights are shown, sum to about 100%, and favour precise studies", {
  p <- plot_forest(eff, feature = "g1")
  lab <- layer(p, "GeomText", 1)$label
  w <- as.numeric(sub("%", "", lab, fixed = TRUE))
  expect_equal(sum(w), 100, tolerance = 0.03)
  se <- c(0.10, 0.20, 0.15, 0.30)                              # Bangladesh, Haiti, USA CA/FL, USA NC
  expect_identical(order(-w), order(se))                       # most weight to the smallest standard error
  sz <- layer(p, "GeomPoint")$size
  expect_equal(max(sz), 4.5)                                   # the largest weight gets the largest marker
})

test_that("studies without an estimate for a feature are labelled", {
  p <- plot_forest(eff, feature = c("g1", "g2"))
  txt <- layer(p, "GeomText", 2)
  expect_equal(nrow(txt), 1)                                   # USA CA/FL has no estimate for g2
  expect_equal(txt$label, "not estimated")
  expect_equal(sum(geoms(plot_forest(eff, feature = "g1")) == "GeomText"), 1)   # nothing missing for g1
})

test_that("a feature seen in too few studies is drawn without a pooled effect", {
  e1 <- bind_effects(mk("s1", "only", 0.3, 0.1), mk("s2", "other", 0.2, 0.1), mk("s3", "other", 0.1, 0.1))
  p <- plot_forest(e1, feature = "only")
  expect_false("GeomPolygon" %in% geoms(p))                                           # no diamond
  expect_match(titles(p), "not pooled", fixed = TRUE)
})

test_that("input is validated", {
  expect_error(plot_forest(data.frame(a = 1), feature = "g1"), "mm_effects")
  expect_error(plot_forest(eff), "feature")
  expect_error(plot_forest(eff, feature = "nope"), "not found in the effect table")
  two <- bind_effects(eff, mk("Bangladesh", "g1", 0.1, 0.1, term = "z"), mk("Haiti", "g1", 0.1, 0.1, term = "z"))
  expect_error(plot_forest(two, feature = "g1"), "choose one with 'term'")
  expect_s3_class(plot_forest(two, feature = "g1", term = "z"), "ggplot")
  expect_error(plot_forest(eff, feature = "g1", level = 2), "level")
  expect_error(plot_forest(eff, pooled = data.frame(a = 1), feature = "g1"), "mm_pooled")
})

test_that("feature labels replace the IDs in the panel titles", {
  ttl <- titles(plot_forest(eff, feature = "g1", feature_labels = c(g1 = "Coprococcus")))
  expect_match(ttl, "^Coprococcus")
})

test_that("the plot can be rendered to a file", {
  skip_on_cran()
  f <- tempfile(fileext = ".png")
  suppressWarnings(ggplot2::ggsave(f, plot_forest(eff, feature = c("g1", "g2")), width = 6, height = 6, dpi = 72))
  expect_true(file.exists(f) && file.size(f) > 1000)
})

# ================ heatmap and volcano plots ================
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
