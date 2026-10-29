# Freezes the outputs of the 1.x function taxa.compare() as reference values for the tests of
# the 2.0 engines (tests/testthat/test-fit-study.R). Needs the 1.x code: run it on the
# `legacy-1.x` branch (or any version that still has taxa.compare() and the dataset taxtab6).
suppressMessages(library(metamicrobiomeR)); load("data-raw/legacy-cran-1.2-data/taxtab6.rda")
tab <- as.data.frame(taxtab6)
phyla <- paste0("k__bacteria.p__", c("actinobacteria", "bacteroidetes", "firmicutes", "proteobacteria"))
tab <- tab[, c("personid", "x.sampleid", "bf", "age.sample", phyla)]
tab$bf <- as.character(tab$bf); tab$personid <- as.character(tab$personid); tab$x.sampleid <- as.character(tab$x.sampleid)
quiet <- function(x) suppressMessages(suppressWarnings(x))
ref <- list(
  tab = tab,
  bezi = quiet(taxa.compare(taxtab = tab, propmed.rel = "gamlss", comvar = "bf", adjustvar = "age.sample", longitudinal = "yes")),
  lmasin_mixed = quiet(taxa.compare(taxtab = tab, propmed.rel = "lm", transform = "asin.sqrt", comvar = "bf", adjustvar = "age.sample", longitudinal = "yes")),
  lmasin_linear = quiet(taxa.compare(taxtab = tab, propmed.rel = "lm", transform = "asin.sqrt", comvar = "bf", adjustvar = "age.sample", longitudinal = "no")),
  note = "Outputs of taxa.compare() from metamicrobiomeR 1.3 (gamlss 5.5-5, lme4) on four phyla of the public Bangladesh data (Subramanian et al. 2014).")
saveRDS(ref, "tests/testthat/fixtures/legacy-reference.rds", version = 2)
cat("saved", format(file.size("tests/testthat/fixtures/legacy-reference.rds") / 1e3, digits = 3), "KB\n")
