# Reproduces the numbers in vignettes/notes-earlier-results.Rmd.
# Args: <results dir> <breastfeeding repo data dir> <dir holding CRAN 1.2 data> [library with legacy_to_study and fit_study]
# results dir: data-raw/legacy-1.1/results as of commit 1e782b2^. Breastfeeding repo: github.com/nhanhocu/metamicrobiome_breastfeeding
a <- commandArgs(TRUE); if (length(a) > 3) .libPaths(c(a[4], .libPaths()))
suppressPackageStartupMessages({library(metamicrobiomeR); library(metafor); library(gamlss)}); options(metamicrobiomeR.quiet = TRUE)
env <- new.env(); for (f in list.files(a[1], full.names = TRUE)) load(f, env)
load(file.path(a[2], "rmdat.shareg7s.rda")); load(file.path(a[3], "gtab.3stud.rda")); load(file.path(a[3], "tabsex4.rda")); load(file.path(a[3], "taxtab6.rda"))
st <- as.data.frame(tabsex4); E <- "Estimate.genderMale"; SE <- "Std. Error.genderMale"
cop <- "k__bacteria.p__firmicutes.c__clostridia.o__clostridiales.f__lachnospiraceae.g__coprococcus"

cat("== Coprococcus ==\n")
u <- as.data.frame(rmdat.unc); u <- u[u$age.sample <= 6, ]; y <- u[[cop]]
cat("NC up to 6 months: samples", nrow(u), "infants", length(unique(u$personid)), "male infants", length(unique(u$personid[u$gender == "Male"])), "\n")
cat("non-zero samples", sum(y > 0), "by sex:", table(u$gender[y > 0]), "\n")
lg <- function(p) log(p / (1 - p)); cat(sprintf("logit difference %.5f, stored estimate %.5f, stored SE %.4f\n", lg(y[y > 0 & u$gender == "Male"]) - lg(y[y > 0 & u$gender == "Female"]), st[st$pop == "USA(NC)" & st$id == cop, E], st[st$pop == "USA(NC)" & st$id == cop, SE]))
pool <- function(d) do.call(rbind, lapply(unique(d$id), function(i) { q <- d[d$id == i & !is.na(d[[E]]) & !is.na(d[[SE]]), ]; if (nrow(q) < 2) return(NULL); f <- suppressWarnings(rma(q[[E]], sei = q[[SE]], method = "DL")); data.frame(id = i, k = nrow(q), est = f$b[1], p = f$pval, maxw = max(weights(f))) }))
sumr <- function(l, d) { P <- pool(d); P$fdr <- p.adjust(P$p, "fdr"); g <- grepl("\\.g__[a-z_\\[\\]-]+$", P$id, perl = TRUE); c0 <- P[P$id == cop, ]; cat(sprintf("%-28s Coprococcus p %.3g, top weight %.1f%% | named genera FDR<.05: %s\n", l, c0$p, c0$maxw, paste(sub(".*g__", "", P$id[g & P$fdr < .05]), collapse = ", "))) }
sumr("all four studies", st); sumr("without USA(NC)", st[st$pop != "USA(NC)", ])

cat("\n== Single-study pooling, stored 1.x result ==\n")
for (tm in names(env$metab.sex$random)) { k <- tapply(!is.na(st[[paste0("Estimate.", tm)]]), st$id, sum); i <- env$metab.sex$random[[tm]]$id; cat(sprintf("%-20s pooled %d, from no study %d, from one study %d\n", tm, length(i), sum(k[i] == 0), sum(k[i] == 1))) }

cat("\n== Non-convergence at 20 iterations, sex model on taxtab6 ==\n")
d <- as.data.frame(unclass(taxtab6), stringsAsFactors = FALSE, check.names = FALSE); tx <- grep("k__", names(d), value = TRUE)
keep <- tx[colSums(d[, tx] > 0, na.rm = TRUE) >= 0.05 * nrow(d) & colMeans(d[, tx], na.rm = TRUE) > 5e-5]
d$personid <- factor(d$personid); d$gender <- factor(d$gender); d$bf <- factor(d$bf)
fit <- function(t, cyc) { td <- d[, c(t, "gender", "bf", "age.sample", "personid")]; names(td)[1] <- "y"; td$y[td$y >= 1] <- .9999; td <- na.omit(td)
  m <- suppressWarnings(gamlss(y ~ gender + bf + age.sample + random(personid), family = BEZI, data = td, trace = FALSE, control = gamlss.control(n.cyc = cyc, trace = FALSE)))
  utils::capture.output(s <- suppressWarnings(summary(m, save = TRUE))); c(conv = m$converged, p = s$coef.table["genderMale", 4]) }
set.seed(20261104); A <- t(sapply(keep, fit, cyc = 20)); nc <- rownames(A)[A[, "conv"] == 0]; cat("not converged:", length(nc), "of", nrow(A), "\n")
sa <- sample(nc, 30); B <- t(sapply(sa, fit, cyc = 300)); cat("sample of 30 at 300 iterations: converged", sum(B[, "conv"]), ", median |change in p|", round(median(abs(A[sa, "p"] - B[, "p"])), 2), ", crossing 0.05:", sum((A[sa, "p"] < .05) != (B[, "p"] < .05)), "\n")

cat("\n== Adjusted p-values: tables against function output, Example 1 ==\n")
s1 <- as.data.frame(env$taxacom6.zi.rmg); ph <- s1[s1$id %in% paste0("k__bacteria.p__", c("actinobacteria", "bacteroidetes", "firmicutes", "proteobacteria")), ]
pc <- "Pr(>|t|).bfNon_exclusiveBF"; ad <- grep("pval.adjust.*Non_exclusiveBF", names(ph), value = TRUE)[1]
print(data.frame(id = sub(".*p__", "", ph$id), p = round(ph[[pc]], 4), stored_function_fdr = round(ph[[ad]], 4), bh_within_4_phyla = round(p.adjust(ph[[pc]], "fdr"), 4)), row.names = FALSE)

cat("\n== Example 3, stored cross-sectional fit against a longitudinal fit ==\n")
for (ex in list(list("EBF <= 2 months", "taxtab6plus.exbf2.rm", "taxacom.6plus.dia.exbf2.zi.rmg"), list("EBF > 2 months", "taxtab6plus.exbf2plus.rm", "taxacom.6plus.dia.exbf2plus.zi.rmg"))) {
  tb <- as.data.frame(get(ex[[2]], env)[[5]]); stored <- as.data.frame(get(ex[[3]], env)); id <- "k__bacteria.p__actinobacteria"
  td <- tb[, c(id, "diarrhea", "age.sample", "personid")]; names(td)[1] <- "y"; td$y[td$y >= 1] <- .9999; td <- na.omit(td); td$personid <- factor(td$personid)
  m <- suppressWarnings(gamlss(y ~ diarrhea + age.sample + random(personid), family = BEZI, data = td, trace = FALSE, control = gamlss.control(n.cyc = 100, trace = FALSE)))
  utils::capture.output(s <- suppressWarnings(summary(m, save = TRUE))); k <- s$coef.table["diarrheaYes", ]
  cat(sprintf("%-16s samples %d, infants %d | stored est %.2f, p %.4f | longitudinal est %.2f, p %.4f\n", ex[[1]], nrow(td), nlevels(td$personid), stored[stored$id == id, "Estimate.diarrheaYes"], stored[stored$id == id, "Pr(>|t|).diarrheaYes"], k[1], k[4]))
}
