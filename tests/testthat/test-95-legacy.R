# The legacy (1.x) functions, checked against stored output of CRAN 1.2 and the stored 2019 pooled result.
# Expected differences are all listed in NEWS.md: the mixed-model p-value factor, the adjusted p-value column names,
# and taxa pooled from too few studies. Fixtures: data-raw/make-legacy-reference.R

ref <- readRDS(test_path("fixtures", "legacy-reference.rds"))
taxtab6 <- metamicrobiomeR::taxtab6
tabsex4 <- as.data.frame(metamicrobiomeR::tabsex4)
tab6 <- as.data.frame(taxtab6)
fuso <- colnames(tab6)[grep("k__bacteria.p__fusobacteria", colnames(tab6))]
sub6 <- tab6[, c("personid", "x.sampleid", "bf", "age.sample", fuso)]
sub_meta <- tabsex4[tabsex4$id %in% tabsex4$id[grep("k__bacteria.p__fusobacteria", tabsex4$id)], ]
tol <- 1e-6
same_cols <- function(new, old) {
  shared <- intersect(colnames(old), colnames(new))
  shared <- shared[vapply(old[shared], is.numeric, logical(1))]
  expect_gt(length(shared), 0)
  expect_equal(as.matrix(new[shared]), as.matrix(old[shared]), tolerance = tol, ignore_attr = TRUE)
}
collect_warnings <- function(expr) {
  w <- character()
  value <- withCallingHandlers(expr, warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  list(value = value, warnings = w)
}

test_that("all legacy functions and data are still there, with their arguments in the same order", {
  for (f in c("alpha.compare", "meta.niceplot", "meta.taxa", "metatab.show", "microbiomeage", "pathway.compare", "read.multi",
              "taxa.compare", "taxa.filter", "taxa.mean.plot", "taxa.meansdn", "taxcomtab.show")) {
    expect_true(is.function(getExportedValue("metamicrobiomeR", f)))
  }
  expect_identical(names(formals(taxa.compare)),
                   c("taxtab", "propmed.rel", "transform", "zeroreplace.method", "comvar", "adjustvar", "personid",
                     "longitudinal", "percent.filter", "relabund.filter", "p.adjust.method", "n.cyc"))
  expect_identical(names(formals(meta.taxa)),
                   c("taxcomdat", "estimate.pattern", "se.pattern", "summary.measure", "pool.var", "studylab",
                     "backtransform", "percent.meta", "p.adjust.method", "method.tau"))
  expect_identical(names(formals(microbiomeage)), c("l6.relabundtab", "bal6", "seed"))
  items <- utils::data(package = "metamicrobiomeR")$results[, "Item"]
  expect_true(all(c("alphadat", "asum4", "covar.rm", "kegg.12", "tabsex4", "taxtab6") %in% items))
  expect_false("gtab.3stud" %in% items)
  expect_false(inherits(taxtab6, "grouped_df"))
})

test_that("taxa.compare() with GAMLSS matches CRAN 1.2 and reports non-convergence once", {
  skip_on_cran()
  r <- collect_warnings(taxa.compare(taxtab = sub6, propmed.rel = "gamlss", comvar = "bf", adjustvar = "age.sample",
                                     longitudinal = "yes", p.adjust.method = "fdr"))
  same_cols(r$value, ref$tc_gamlss_long)
  expect_identical(r$value$id, ref$tc_gamlss_long$id)
  info <- attr(r$value, "fit_info")
  expect_named(info, c("id", "n", "n_nonzero", "n_coef", "converged"))
  expect_equal(sum(grepl("did not converge within n.cyc = 20", r$warnings)), as.integer(any(!info$converged)))
  expect_false(any(grepl("has not yet converged", r$warnings)))
  cs <- suppressWarnings(taxa.compare(taxtab = sub6, propmed.rel = "gamlss", comvar = "bf", adjustvar = "age.sample",
                                      longitudinal = "no", p.adjust.method = "fdr"))
  same_cols(cs, ref$tc_gamlss_cs)
})

test_that("n.cyc changes the GAMLSS iterations", {
  skip_on_cran()
  one <- sub6[, c("personid", "bf", "age.sample", "k__bacteria.p__fusobacteria")]
  r <- suppressWarnings(taxa.compare(one, propmed.rel = "gamlss", comvar = "bf", adjustvar = "age.sample", longitudinal = "yes", n.cyc = 300))
  expect_true(attr(r, "fit_info")$converged)
})

test_that("mixed-model p-values are two-sided, 1.2 values were 1.96/2 of them, and nothing else changes", {
  new <- suppressMessages(taxa.compare(taxtab = sub6, propmed.rel = "lm", transform = "asin.sqrt", comvar = "bf", adjustvar = "age.sample", longitudinal = "yes"))
  old <- ref$tc_lm_long
  est <- grep("^(Estimate|Std\\. Error|t value)", colnames(old), value = TRUE)
  expect_equal(as.matrix(new[est]), as.matrix(old[est]), tolerance = tol, ignore_attr = TRUE)
  p <- grep("^Pr", colnames(old), value = TRUE)
  expect_equal(as.numeric(as.matrix(old[p]) / as.matrix(new[p])), rep(1.96 / 2, nrow(old) * length(p)), tolerance = tol)
  cs <- suppressMessages(taxa.compare(taxtab = sub6, propmed.rel = "lm", transform = "asin.sqrt", comvar = "bf", adjustvar = "age.sample", longitudinal = "no"))
  same_cols(cs, ref$tc_lm_cs)
})

test_that("adjusted p-values have full names and the old truncated name stays as an alias", {
  new <- suppressMessages(taxa.compare(taxtab = sub6, propmed.rel = "lm", transform = "asin.sqrt", comvar = "bf", adjustvar = "age.sample", longitudinal = "no"))
  expect_true(all(c("pval.adjust.age.sample", "pval.adjust.sample") %in% colnames(new)))
  expect_equal(new$pval.adjust.age.sample, new$pval.adjust.sample)
})

test_that("taxa with too few non-zero samples are flagged", {
  set.seed(1)
  d <- data.frame(personid = rep(paste0("p", 1:6), each = 3), group = rep(c("a", "b"), each = 9), age.sample = stats::runif(18, 0, 6))
  d$k__bacteria.p__rare <- 0
  d$k__bacteria.p__rare[c(2, 11)] <- c(0.0008, 0.0012)
  d$k__bacteria.p__common <- stats::rbeta(18, 2, 20)
  r <- collect_warnings(taxa.compare(d, propmed.rel = "gamlss", comvar = "group", adjustvar = "age.sample", longitudinal = "no"))
  info <- attr(r$value, "fit_info")
  expect_equal(info$n_nonzero[info$id == "k__bacteria.p__rare"], 2)
  expect_equal(info$n_coef[info$id == "k__bacteria.p__rare"], 3)
  expect_true(any(grepl("no more non-zero samples than mean-model coefficients.*k__bacteria.p__rare", r$warnings)))
  expect_false(any(grepl("p__common", r$warnings)))
})

test_that("the logit transform leaves out taxa with zeros and says so", {
  set.seed(2)
  d <- data.frame(personid = paste0("p", 1:12), group = rep(c("a", "b"), 6), age.sample = 1:12)
  d$k__bacteria.p__zero <- c(0, stats::runif(11, 0.01, 0.2))
  d$k__bacteria.p__nozero <- stats::runif(12, 0.01, 0.2)
  expect_warning(new <- taxa.compare(d, propmed.rel = "lm", transform = "logit", comvar = "group", adjustvar = "age.sample", longitudinal = "no"), "1 of 2 taxa have zeros")
  expect_identical(new$id, "k__bacteria.p__nozero")
  expect_error(taxa.compare(d[, c("personid", "group", "age.sample", "k__bacteria.p__zero")], propmed.rel = "lm", transform = "logit",
                            comvar = "group", adjustvar = "age.sample", longitudinal = "no"), "every taxon has zeros")
})

test_that("the GMPR transform gives usable results", {
  tl <- colnames(tab6)[grep("^k__bacteria\\.p__[a-z]+$", colnames(tab6))]
  new <- suppressWarnings(suppressMessages(taxa.compare(tab6[, c("personid", "bf", "age.sample", tl)], propmed.rel = "lm", transform = "gmpr",
                                                        comvar = "bf", adjustvar = "age.sample", longitudinal = "no")))
  expect_gt(nrow(new), 0)
  expect_true(all(is.finite(new[["Estimate.bfNon_exclusiveBF"]])))
})

test_that("meta.taxa() reproduces the stored 2019 pooled result on the full table", {
  new <- suppressWarnings(meta.taxa(taxcomdat = tabsex4, summary.measure = "RR", pool.var = "id", studylab = "study",
                                    backtransform = FALSE, percent.meta = 0.5, p.adjust.method = "fdr"))
  for (w in c("random", "fixed")) {
    o <- ref$stored_2019_sex[[w]]
    n <- new[[w]]$genderMale
    expect_equal(nrow(n), 185)
    expect_setequal(n$id, o$id)
    m <- merge(o, n, by = "id", suffixes = c(".o", ".n"))
    for (v in c("estimate", "se", "ll", "ul", "statistic", "p")) expect_equal(m[[paste0(v, ".n")]], m[[paste0(v, ".o")]], tolerance = 1e-8)
  }
})

test_that("meta.taxa() needs two studies with an estimate, flags dominance and accepts another estimator", {
  r <- collect_warnings(meta.taxa(taxcomdat = tabsex4, summary.measure = "RR", pool.var = "id", studylab = "study", percent.meta = 0.5))
  new <- r$value
  k <- tapply(!is.na(tabsex4[["Estimate.bfNo_BF"]]), tabsex4$id, sum)
  expect_true(all(k[new$random$bfNo_BF$id] >= 2))
  expect_equal(nrow(new$random$bfNo_BF), sum(k >= 2))
  expect_true(any(grepl("bfNo_BF: .* taxa were not pooled", r$warnings)))
  dom <- attr(new, "dominance")
  cop <- dom[dom$term == "genderMale" & grepl("g__coprococcus$", dom$id), ]
  expect_match(cop$top_study, "USA(NC)", fixed = TRUE)
  expect_gt(cop$max_weight, 0.99)
  reml <- suppressWarnings(meta.taxa(sub_meta, method.tau = "REML"))
  expect_equal(colnames(reml$random$genderMale), c("estimate", "se", "ll", "ul", "statistic", "p", "id", "p.adjust"))
})

test_that("the other legacy functions match CRAN 1.2", {
  skip_on_cran()
  kegg.12 <- metamicrobiomeR::kegg.12; covar.rm <- metamicrobiomeR::covar.rm; alphadat <- metamicrobiomeR::alphadat
  expect_identical(taxa.filter(taxtab = tab6, percent.filter = 0.05, relabund.filter = 0.00005), ref$filter)
  same_cols(as.data.frame(taxa.meansdn(taxtab = tab6, sumvar = "bf", groupvar = "age.sample")), ref$meansdn)
  p <- suppressWarnings(pathway.compare(pathtab = list(kegg.12[[1]][, 1:2]), mapfile = covar.rm, sampleid = "sampleid", pathsum = "rel", stat.med = "gamlss",
                                        comvar = "gender", adjustvar = c("age.sample", "bf"), longitudinal = "no", p.adjust.method = "fdr",
                                        percent.filter = 0.05, relabund.filter = 0.00005))$l1
  same_cols(p, ref$path)
  cv <- covar.rm
  cv$sampleid <- tolower(cv$sampleid)
  a <- suppressWarnings(suppressMessages(alpha.compare(datlist = alphadat, depth = 3, mapfile = cv, mapsampleid = "sampleid", comvar = "gender",
                                                       adjustvar = c("age.sample", "bf"), longitudinal = "yes", age.limit = 6, standardize = TRUE)))$alphasum
  same_cols(a, ref$alpha)
  m <- suppressWarnings(meta.taxa(taxcomdat = sub_meta, summary.measure = "RR", pool.var = "id", studylab = "study", backtransform = FALSE,
                                  percent.meta = 0.5, p.adjust.method = "fdr"))
  same_cols(m$random$genderMale, ref$meta_sub_random)
  same_cols(m$fixed$genderMale, ref$meta_sub_fixed)
  tab <- metatab.show(metatab = m$random, com.pooled.tab = sub_meta, tax.lev = "l2", showvar = "genderMale", p.cutoff.type = "p", p.cutoff = 1, display = "table")
  same_cols(tab, ref$metatab)
  same_cols(taxcomtab.show(taxcomtab = p, sumvar = "path", tax.lev = "l2", tax.select = "none", showvar = "genderMale", p.adjust.method = "fdr", p.cutoff = 1), ref$taxcomtab)
})

test_that("read.multi() reads a folder", {
  d <- tempfile()
  dir.create(d)
  utils::write.table(data.frame(A = 1), file.path(d, "one.txt"), sep = "\t", row.names = FALSE)
  r <- read.multi(patht = d, patternt = ".txt", study = "s")
  expect_named(r, "one")
  expect_equal(r$one$study, "s")
  expect_named(r$one, c("a", "study"))
})

test_that("legacy_to_study() moves a table to the modern workflow and gives the same GAMLSS estimate", {
  skip_on_cran()
  s <- suppressWarnings(legacy_to_study(sub6, list(id = "Bangladesh", feature_type = "otu", subject = "personid")))
  expect_s3_class(s, "mm_study")
  expect_equal(mm_features(s)["k__bacteria.p__fusobacteria", "rank"], "phylum")
  f <- suppressWarnings(suppressMessages(fit_study(s, ~ bf + age.sample, rank = "phylum", longitudinal = TRUE, n_cycles = 20, exclude_unassigned = FALSE)))
  f <- f[f$term == "bfNon_exclusiveBF", ]
  old <- ref$tc_gamlss_long
  old <- old[old$id == "k__bacteria.p__fusobacteria", ]
  expect_equal(f$estimate, old[["Estimate.bfNon_exclusiveBF"]], tolerance = 1e-8)
  expect_error(legacy_to_study(sub6, list(id = "x", feature_type = "otu"), id_col = "nope"), "is not in the table")
})
