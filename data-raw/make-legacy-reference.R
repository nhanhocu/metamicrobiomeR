# Generates tests/testthat/fixtures/legacy-reference.rds from an installed CRAN metamicrobiomeR 1.2
# and the stored 2019 pooled result. Args: <library holding 1.2> <metab.sex.rda> <output rds>
a <- commandArgs(TRUE); .libPaths(c(a[1], .libPaths())); suppressPackageStartupMessages({library(metamicrobiomeR); library(meta)})
stopifnot(packageVersion("metamicrobiomeR") == "1.2")
q <- function(expr) { utils::capture.output(v <- suppressWarnings(expr)); v }
data(taxtab6); tab6 <- as.data.frame(taxtab6); tl <- colnames(taxtab6)[grep("k__bacteria.p__fusobacteria", colnames(taxtab6))]
sub6 <- tab6[, c("personid", "x.sampleid", "bf", "age.sample", tl)]
r <- list()
r$tc_gamlss_long <- q(taxa.compare(taxtab = sub6, propmed.rel = "gamlss", comvar = "bf", adjustvar = "age.sample", longitudinal = "yes", p.adjust.method = "fdr"))
r$tc_gamlss_cs <- q(taxa.compare(taxtab = sub6, propmed.rel = "gamlss", comvar = "bf", adjustvar = "age.sample", longitudinal = "no", p.adjust.method = "fdr"))
r$tc_lm_long <- q(taxa.compare(taxtab = sub6, propmed.rel = "lm", transform = "asin.sqrt", comvar = "bf", adjustvar = "age.sample", longitudinal = "yes"))
r$tc_lm_cs <- q(taxa.compare(taxtab = sub6, propmed.rel = "lm", transform = "asin.sqrt", comvar = "bf", adjustvar = "age.sample", longitudinal = "no"))
r$filter <- taxa.filter(taxtab = tab6, percent.filter = 0.05, relabund.filter = 0.00005)
r$meansdn <- as.data.frame(q(taxa.meansdn(taxtab = tab6, sumvar = "bf", groupvar = "age.sample")))
data(kegg.12); data(covar.rm)
r$path <- q(pathway.compare(pathtab = list(kegg.12[[1]][, 1:2]), mapfile = covar.rm, sampleid = "sampleid", pathsum = "rel", stat.med = "gamlss", comvar = "gender", adjustvar = c("age.sample", "bf"), longitudinal = "no", p.adjust.method = "fdr", percent.filter = 0.05, relabund.filter = 0.00005))$l1
data(alphadat); cv <- covar.rm; cv$sampleid <- tolower(cv$sampleid)
r$alpha <- q(alpha.compare(datlist = alphadat, depth = 3, mapfile = cv, mapsampleid = "sampleid", comvar = "gender", adjustvar = c("age.sample", "bf"), longitudinal = "yes", age.limit = 6, standardize = TRUE))$alphasum
data(tabsex4); tlm <- tabsex4$id[grep("k__bacteria.p__fusobacteria", tabsex4$id)]
settings.meta(method.tau = "DL")
m <- q(meta.taxa(taxcomdat = tabsex4[tabsex4$id %in% tlm, ], summary.measure = "RR", pool.var = "id", studylab = "study", backtransform = FALSE, percent.meta = 0.5, p.adjust.method = "fdr"))
r$meta_sub_random <- m$random$genderMale; r$meta_sub_fixed <- m$fixed$genderMale
r$metatab <- metatab.show(metatab = m$random, com.pooled.tab = tabsex4[tabsex4$id %in% tlm, ], tax.lev = "l2", showvar = "genderMale", p.cutoff.type = "p", p.cutoff = 1, display = "table")
r$taxcomtab <- q(taxcomtab.show(taxcomtab = r$path, sumvar = "path", tax.lev = "l2", tax.select = "none", showvar = "genderMale", p.adjust.method = "fdr", p.cutoff = 1))
e <- new.env(); load(a[2], e)
r$stored_2019_sex <- list(random = e$metab.sex$random$genderMale, fixed = e$metab.sex$fixed$genderMale)
r$versions <- list(R = R.version.string, gamlss = as.character(packageVersion("gamlss")), meta = as.character(packageVersion("meta")), lme4 = as.character(packageVersion("lme4")))
saveRDS(r, a[3], version = 2, compress = "xz"); cat("saved", length(r), "objects\n"); print(sapply(r, function(x) if (is.data.frame(x)) paste(dim(x), collapse = "x") else class(x)[1]))
