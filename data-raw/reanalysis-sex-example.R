# Re-analysis of the 1.x sex-difference example with the functions of version 2.0, for the paper.
#
# Usage: Rscript data-raw/reanalysis-sex-example.R <breastfeeding-data-dir> <cran-1.2-data-dir> <output-dir> [cores]
#   <breastfeeding-data-dir>: the data folder of https://github.com/nhanhocu/metamicrobiome_breastfeeding (rmdat.shareg7s.rda)
#   <cran-1.2-data-dir>: the data folder of metamicrobiomeR 1.2 from the CRAN archive (gtab.3stud.rda, taxtab6.rda)
#
# Infants up to 6 months in four cohorts. Per cohort: relative abundance ~ gender + bf + age.sample, with a subject random
# intercept where infants have repeated samples (not Haiti). The estimates for genderMale are pooled across cohorts.
# Runs: two engines (bezi_gamlss, bezi_glmmtmb), each with the default low-support guard and with min_nonzero = 0.

suppressPackageStartupMessages(library(metamicrobiomeR))
a <- commandArgs(TRUE)
bfd <- a[1]; crd <- a[2]; out <- a[3]; options(metamicrobiomeR.cores = if (length(a) >= 4) as.integer(a[4]) else 2L)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

e <- new.env()
load(file.path(bfd, "rmdat.shareg7s.rda"), e); load(file.path(crd, "gtab.3stud.rda"), e); load(file.path(crd, "taxtab6.rda"), e)

# Cohort tables in the 1.x layout: mapping variables and k__ taxa in columns, one row per sample
mk <- function(rm, g, prefix, pid) {
  u <- as.data.frame(get(rm, e)); gg <- as.data.frame(e$gtab.3stud[[g]])
  u$key <- sub(prefix, "", tolower(u$sampleid)); gg$key <- tolower(as.character(gg$x.sampleid))
  m <- u[, c("key", if (pid) "personid", "gender", "bf", "age.sample")]
  d <- merge(m, gg[, setdiff(names(gg), c("age.sample", "x.sampleid"))], by = "key")
  d <- d[d$age.sample <= 6 & !is.na(d$gender) & !is.na(d$bf), ]
  d$x.sampleid <- d$key; d$key <- NULL
  d
}
cohorts <- list(
  `USA(NC)` = list(tab = mk("rmdat.unc", "nc", "^unc\\.unc\\.", TRUE), long = TRUE),
  `USA(CA_FL)` = list(tab = mk("rmdat.usbmk", "ca_fl", "^usbmk\\.usbmk\\.", TRUE), long = TRUE),
  Haiti = list(tab = mk("rmdat.ha", "haiti", "^ha\\.ha\\.", FALSE), long = FALSE),
  Bangladesh = list(tab = as.data.frame(unclass(e$taxtab6), stringsAsFactors = FALSE, check.names = FALSE), long = TRUE)
)
studies <- lapply(names(cohorts), function(n) {
  x <- cohorts[[n]]$tab
  x$gender <- factor(as.character(x$gender), levels = c("Female", "Male"))
  x$bf <- droplevels(factor(x$bf))
  meta <- list(id = n, feature_type = "otu", subject = if (cohorts[[n]]$long) "personid")
  suppressMessages(legacy_to_study(x, meta, id_col = "x.sampleid"))
})
names(studies) <- names(cohorts)
desc <- do.call(rbind, lapply(names(studies), function(n) {
  s <- mm_samples(studies[[n]])
  data.frame(cohort = n, samples = nrow(s), infants = if (cohorts[[n]]$long) length(unique(s$personid)) else nrow(s),
             female = sum(s$gender == "Female"), male = sum(s$gender == "Male"), features = nrow(mm_relabund(studies[[n]])))
}))
print(desc)
utils::write.csv(desc, file.path(out, "cohorts.csv"), row.names = FALSE)

runs <- expand.grid(engine = c("bezi_gamlss", "bezi_glmmtmb"), guard = c(TRUE, FALSE), stringsAsFactors = FALSE)
pooled <- list(); fitted <- list()
for (i in seq_len(nrow(runs))) {
  r <- runs[i, ]; lab <- paste0(r$engine, if (r$guard) "_guard" else "_noguard")
  fits <- lapply(names(studies), function(n) {
    f <- suppressWarnings(suppressMessages(fit_study(studies[[n]], ~ gender + bf + age.sample, engine = r$engine,
                                                     longitudinal = cohorts[[n]]$long, min_prevalence = 0.05,
                                                     min_nonzero = if (r$guard) NULL else 0)))
    ls <- attr(f, "low_support")
    cat(lab, n, ": fitted", length(unique(f$feature)), "| left out for low support", if (is.null(ls)) 0 else nrow(ls), "\n")
    f
  })
  eff <- suppressMessages(do.call(bind_effects, fits))
  fitted[[lab]] <- eff
  for (pm in c("REML_adhoc", "DL_z")) {
    p <- suppressWarnings(suppressMessages(pool_effects(eff, term = "genderMale", part = "abundance",
                                                         method = if (pm == "DL_z") "DL" else "REML",
                                                         test = if (pm == "DL_z") "z" else "adhoc")))
    p$run <- lab; p$pooling <- pm
    pooled[[paste(lab, pm)]] <- p
  }
}
P <- do.call(rbind, lapply(pooled, as.data.frame))
utils::write.csv(P, file.path(out, "pooled.csv"), row.names = FALSE)
saveRDS(list(fitted = fitted, pooled = P, cohorts = desc), file.path(out, "reanalysis.rds"))

# Summary: Coprococcus and the taxa significant after FDR correction, for each run and pooling method
cop <- grepl("g__coprococcus$", P$feature)
summ <- do.call(rbind, lapply(split(P, list(P$run, P$pooling), drop = TRUE), function(d) {
  cp <- d[grepl("g__coprococcus$", d$feature), ]
  sig <- d$feature[which(d$p_adjust < 0.05)]
  data.frame(run = d$run[1], pooling = d$pooling[1], taxa_pooled = nrow(d),
             coprococcus_k = if (nrow(cp)) cp$k else NA, coprococcus_p = if (nrow(cp)) signif(cp$p_value, 3) else NA,
             coprococcus_max_weight = if (nrow(cp)) round(cp$max_weight, 3) else NA,
             significant_fdr = length(sig), which = paste(sub(".*\\.", "", sig), collapse = "; "))
}))
print(summ, row.names = FALSE)
utils::write.csv(summ, file.path(out, "summary.csv"), row.names = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "sessionInfo.txt"))
