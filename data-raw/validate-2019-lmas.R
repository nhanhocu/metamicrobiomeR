# Reproduce the LMAS (linear model on arcsine square root abundance) columns of Tables 3-5
# of Ho et al. (2019) with the 2.0 `lm_asin` engine, against the stored 1.x results and the
# printed numbers. Run from the repository root: Rscript data-raw/validate-2019-lmas.R
source("data-raw/legacy-helpers.R")
env <- new.env(); for (f in list.files("data-raw/legacy-1.1/results", full.names = TRUE)) load(f, env)
phyla <- c("actinobacteria", "bacteroidetes", "firmicutes", "proteobacteria")
P <- function(ex, est, lo, hi, p, fdr) data.frame(example = ex, phylum = phyla, est, lo, hi, p, fdr)
paper <- rbind(
  P("1",  c(-0.13, 0.03, 0.07, 0.10), c(-0.23, 0.00, 0.00, 0.02), c(-0.03, 0.05, 0.14, 0.17), c(.0088, .0292, .0668, .0103), c(.0207, .0390, .0668, .0207)),
  P("2",  c(0.05, -0.05, -0.04, 0.02), c(-0.06, -0.09, -0.12, -0.02), c(0.16, -0.01, 0.04, 0.07), c(.3451, .0270, .3168, .2916), c(.3451, .1079, .3451, .3451)),
  P("3a", c(-0.12, 0.06, 0.11, 0.00), c(-0.23, -0.12, 0.01, -0.07), c(0.00, 0.01, 0.20, 0.08), c(.0424, .0852, .0269, .9060), c(.0848, .1136, .0848, .9060)),
  P("3b", c(0.00, 0.01, -0.01, 0.02), c(-0.10, -0.07, -0.13, -0.06), c(0.10, 0.09, 0.12, 0.11), c(.9626, .8101, .8927, .5875), c(.9989, .9707, .9707, .9191)))
examples <- list(
  "1"  = list(tab = "taxtab6.rm",            f = ~ bf + age.sample,          term = "bfNon_exclusiveBF",      long = TRUE,  stored = "taxacom6.rmg.as"),
  "2"  = list(tab = "taxtab6plus.rm",        f = ~ month.food5 + age.sample, term = "month.food5>5 months",   long = TRUE,  stored = "taxacom.6plus.sl5.rmg.as"),
  "3a" = list(tab = "taxtab6plus.exbf2.rm",     f = ~ diarrhea + age.sample, term = "diarrheaYes",            long = FALSE, stored = "taxacom.6plus.dia.exbf2.rmg.as"),
  "3b" = list(tab = "taxtab6plus.exbf2plus.rm", f = ~ diarrhea + age.sample, term = "diarrheaYes",            long = FALSE, stored = "taxacom.6plus.dia.exbf2plus.rmg.as"))
short <- function(x) sub("^k__bacteria\\.p__", "", x)
for (ex in names(examples)) {
  e <- examples[[ex]]
  s <- legacy_to_study(as.data.frame(get(e$tab, env)[[5]]), list(id = ex, feature_type = "otu", subject = "personid"))$study
  fit <- fit_study(s, e$f, engine = "lm_asin", rank = "phylum", longitudinal = e$long, exclude_unassigned = FALSE)
  fit <- fit[fit$term == e$term & short(fit$feature) %in% phyla, ]; fit$phylum <- short(fit$feature)
  fit$lo <- fit$estimate - 1.96 * fit$se; fit$hi <- fit$estimate + 1.96 * fit$se
  fit$fdr <- p.adjust(fit$p_value, "fdr")
  st <- get(e$stored, env); st$phylum <- short(st$id); m <- merge(fit, st, by = "phylum")
  ratio <- m$p_value / m[[paste0("Pr(>|t|).", e$term)]]
  cat(sprintf("Example %-2s (%s): max |diff| vs stored: estimate %.1e, se %.1e | new p / stored p: %s\n", ex, if (e$long) "mixed model" else "linear model",
      max(abs(m$estimate - m[[paste0("Estimate.", e$term)]])), max(abs(m$se - m[[paste0("Std. Error.", e$term)]])), paste(round(ratio, 4), collapse = ", ")))
  pp <- paper[paper$example == ex, ]; names(pp)[3:7] <- paste0("paper_", names(pp)[3:7])
  nw <- fit[, c("phylum", "estimate", "lo", "hi", "p_value", "fdr")]; names(nw)[2:6] <- c("new_est", "new_lo", "new_hi", "new_p", "new_fdr")
  g <- merge(pp, nw, by = "phylum")
  f98 <- if (e$long) 0.98 else 1                       # the paper's mixed-model p-values carry the 1.96*pnorm error
  g$p_as_published <- g$new_p * f98
  ok <- cbind(est = abs(round(g$new_est, 2) - g$paper_est) < 0.0051, lo = abs(round(g$new_lo, 2) - g$paper_lo) < 0.0051,
              hi = abs(round(g$new_hi, 2) - g$paper_hi) < 0.0051, p = abs(round(g$p_as_published, 4) - g$paper_p) < 0.00005)
  cat(sprintf("           vs paper (p scaled by %.2f): estimate %d/4, lower %d/4, upper %d/4, p %d/4\n", f98, sum(ok[, "est"]), sum(ok[, "lo"]), sum(ok[, "hi"]), sum(ok[, "p"])))
  if (!all(ok)) print(data.frame(phylum = g$phylum, new_est = round(g$new_est, 3), paper_est = g$paper_est, new_lo = round(g$new_lo, 3), paper_lo = g$paper_lo,
                                 new_hi = round(g$new_hi, 3), paper_hi = g$paper_hi, p_pub = round(g$p_as_published, 4), paper_p = g$paper_p)[!apply(ok, 1, all), ], row.names = FALSE)
  g$fdr_pub <- p.adjust(g$p_as_published, "fdr")      # BH over the four phyla, from the p-values as published
  fd <- abs(round(g$fdr_pub, 4) - g$paper_fdr) < 0.00005
  g$new_fdr <- g$fdr_pub
  cat(sprintf("           FDR (BH within the 4 phyla, from the p-values as published): %d/4 match the paper%s\n", sum(fd), if (all(fd)) "" else paste0("; differing: ", paste(g$phylum[!fd], round(g$new_fdr[!fd], 4), "vs", g$paper_fdr[!fd], collapse = "; "))))
}

# Example 3b: are two phylum rows of the paper's table swapped?
e <- examples[["3b"]]
s <- legacy_to_study(as.data.frame(get(e$tab, env)[[5]]), list(id = "3b", feature_type = "otu", subject = "personid"))$study
fit <- fit_study(s, e$f, engine = "lm_asin", rank = "phylum", longitudinal = FALSE, exclude_unassigned = FALSE)
fit <- fit[fit$term == e$term, ]; fit$phylum <- short(fit$feature)
pp <- paper[paper$example == "3b", ]
cat("\n== Example 3b, computed rows vs the paper's rows (matched by p-value) ==\n")
cmp <- merge(data.frame(phylum_paper = pp$phylum, paper_est = pp$est, paper_p = pp$p), fit[fit$phylum %in% phyla, c("phylum", "estimate", "p_value")], by = NULL)
cmp <- cmp[abs(cmp$paper_p - round(cmp$p_value, 4)) < 0.00006, ]
print(data.frame(paper_label = cmp$phylum_paper, computed_label = cmp$phylum, paper_est = cmp$paper_est, computed_est = round(cmp$estimate, 3), p = round(cmp$p_value, 4)), row.names = FALSE)
allp <- fit$p_value; cat("BH over all", length(allp), "phylum-level features, for the four displayed phyla:", paste(round(p.adjust(allp, "fdr")[fit$phylum %in% phyla], 4), collapse = ", "), "| paper:", paste(pp$fdr, collapse = ", "), "\n")
