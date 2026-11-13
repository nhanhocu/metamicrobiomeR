# Fixes to the legacy (1.x) functions for problems reported in GitHub issues. Each test names the issue.

tab6 <- as.data.frame(metamicrobiomeR::taxtab6)
phyla <- c("k__bacteria.p__actinobacteria", "k__bacteria.p__bacteroidetes", "k__bacteria.p__firmicutes", "k__bacteria.p__proteobacteria")
d <- tab6[, c("personid", "gender", "age.sample", "bf", phyla)]
d$chimio <- factor(ifelse(d$gender == "Male", "Pre", "Post"), levels = c("Post", "Pre"))
d$grp <- factor(ifelse(d$gender == "Male", "B", "A"))

test_that("a covariate level containing 'Pr' does not change the estimates (#18)", {
  skip_if_not_installed("gamlss")
  a <- suppressWarnings(taxa.compare(taxtab = d, propmed.rel = "gamlss", comvar = "chimio", adjustvar = "age.sample", longitudinal = "no"))
  b <- suppressWarnings(taxa.compare(taxtab = d, propmed.rel = "gamlss", comvar = "grp", adjustvar = "age.sample", longitudinal = "no"))
  expect_equal(a$Estimate.chimioPre, b$Estimate.grpB, tolerance = 1e-8)
  expect_equal(a[["Std. Error.chimioPre"]], b[["Std. Error.grpB"]], tolerance = 1e-8)
  expect_equal(a$pval.adjust.chimioPre, b$pval.adjust.grpB, tolerance = 1e-8)
  expect_false(any(grepl("\\.[0-9]+$", colnames(a))))
  expect_true(all(a$pval.adjust.chimioPre >= 0 & a$pval.adjust.chimioPre <= 1))
})

test_that("taxa.compare() and pathway.compare() run without adjustment variables (#17)", {
  r <- suppressWarnings(taxa.compare(taxtab = d, propmed.rel = "lm", comvar = "gender", longitudinal = "no"))
  expect_s3_class(r, "data.frame")
  expect_true("Estimate.genderMale" %in% colnames(r))
  expect_null(formals(pathway.compare)$adjustvar)
})

test_that("taxon names that are not syntactic R names are fitted (#5)", {
  skip_if_not_installed("gamlss")
  d5 <- d
  colnames(d5)[colnames(d5) == phyla[2]] <- "k__Bacteria; p__Bacteroidetes-1"
  for (m in c("lm", "gamlss")) {
    r <- suppressWarnings(taxa.compare(taxtab = d5, propmed.rel = m, comvar = "gender", adjustvar = "age.sample", longitudinal = "no"))
    expect_true("k__Bacteria; p__Bacteroidetes-1" %in% r$id, info = m)
  }
})

test_that("names starting with d__ give a clear message (#14)", {
  d14 <- d
  colnames(d14) <- sub("^k__", "d__", colnames(d14))
  expect_error(taxa.compare(taxtab = d14, propmed.rel = "lm", comvar = "gender", longitudinal = "no"), "rename them first")
  expect_error(taxa.filter(d14), "No taxon columns found")
})

test_that("tax.select accepts several taxa, and a table with one taxonomic level is shown (#7)", {
  b <- suppressWarnings(taxa.compare(taxtab = d, propmed.rel = "lm", comvar = "grp", adjustvar = "age.sample", longitudinal = "no"))
  t2 <- taxcomtab.show(taxcomtab = b, tax.lev = "l2", tax.select = phyla[1:2], showvar = "grpB", p.cutoff = 1)
  expect_setequal(t2$id, phyla[1:2])
  expect_equal(nrow(taxcomtab.show(taxcomtab = b, tax.lev = "l2", showvar = "grpB", p.cutoff = 1)), length(phyla))
})

test_that("taxa.mean.plot() accepts several taxa in tax.select and shows phyla only at l2 (#3)", {
  m <- taxa.meansdn(taxtab = tab6, sumvar = "bf", groupvar = "age.sample")
  p <- taxa.mean.plot(tabmean = m, tax.lev = "l2", tax.select = phyla[1:3], comvar = "bf", groupvar = "age.sample")
  expect_setequal(p$taxuse.rm, phyla[1:3])
  m$`k__archaea_mean` <- 0.01
  p2 <- taxa.mean.plot(tabmean = m, tax.lev = "l2", comvar = "bf", groupvar = "age.sample", mean.filter = 0.005)
  expect_false("k__archaea" %in% p2$taxuse.rm)
  expect_true(all(grepl("p__", p2$taxuse.rm)))
})

test_that("taxa.meansdn() accepts several variables in othervar and leaves out sumvar and groupvar (#2)", {
  dd <- d[, c("gender", "age.sample", "bf", phyla)]
  r <- taxa.meansdn(taxtab = dd, sumvar = "bf", groupvar = "age.sample", othervar = c("gender"))
  expect_false(any(grepl("^(sumvar|groupvar|bf|age.sample)_mean$", colnames(r))))
  expect_true(all(paste0(phyla, "_mean") %in% colnames(r)))
})

test_that("metatab.show() works without a 'pop' column and explains an empty result (#19, #22)", {
  tabsex4 <- as.data.frame(metamicrobiomeR::tabsex4)
  sub_meta <- tabsex4[tabsex4$id %in% tabsex4$id[grep("k__bacteria.p__fusobacteria", tabsex4$id)], ]
  m <- suppressWarnings(meta.taxa(taxcomdat = sub_meta, summary.measure = "RR", pool.var = "id", studylab = "study", backtransform = FALSE,
                                  percent.meta = 0.5, p.adjust.method = "fdr"))
  nopop <- sub_meta[, colnames(sub_meta) != "pop"]
  a <- metatab.show(metatab = m$random, com.pooled.tab = sub_meta, tax.lev = "l2", showvar = "genderMale", p.cutoff = 1, display = "data")
  b <- metatab.show(metatab = m$random, com.pooled.tab = nopop, tax.lev = "l2", showvar = "genderMale", p.cutoff = 1, display = "data")
  expect_equal(nrow(b$taxsig.all), nrow(a$taxsig.all))
  expect_equal(b$taxsig.all$estimate, a$taxsig.all$estimate)
  expect_warning(metatab.show(metatab = m$random, com.pooled.tab = sub_meta, tax.lev = "l7", highest.lev = "s", showvar = "genderMale",
                              p.cutoff = 1, display = "table"), "No pooled result at tax.lev")
  expect_no_warning(metatab.show(metatab = m$random, com.pooled.tab = sub_meta, tax.lev = "l2", showvar = "genderMale", p.cutoff = 0,
                                 display = "table"))
})

test_that("GAMLSS without a subject effect does not need a personid column (#9)", {
  skip_if_not_installed("gamlss")
  d9 <- d[, colnames(d) != "personid"]
  r <- suppressWarnings(taxa.compare(taxtab = d9, propmed.rel = "gamlss", comvar = "gender", adjustvar = c("age.sample", "bf"), longitudinal = "no"))
  expect_s3_class(r, "data.frame")
  expect_true("Estimate.genderMale" %in% colnames(r))
})
