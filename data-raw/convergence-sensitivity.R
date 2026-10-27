# Do the 2019 results depend on GAMLSS fits that hit the default iteration limit (20 cycles)?
# Refits the four main phyla of each example with a much higher limit and compares.
source("data-raw/legacy-helpers.R")
env <- new.env(); for (f in list.files("data-raw/legacy-1.1/results", full.names = TRUE)) load(f, env)
phyla <- paste0("k__bacteria.p__", c("actinobacteria", "bacteroidetes", "firmicutes", "proteobacteria"))
ex <- list(
  "1"  = list(tab = "taxtab6.rm",            f = ~ bf + age.sample,          term = "bfNon_exclusiveBF",      long = TRUE),
  "2"  = list(tab = "taxtab6plus.rm",        f = ~ month.food5 + age.sample, term = "month.food5>5 months",   long = TRUE),
  "3a" = list(tab = "taxtab6plus.exbf2.rm",     f = ~ diarrhea + age.sample, term = "diarrheaYes",            long = FALSE),
  "3b" = list(tab = "taxtab6plus.exbf2plus.rm", f = ~ diarrhea + age.sample, term = "diarrheaYes",            long = FALSE))
out <- list()
for (k in names(ex)) {
  e <- ex[[k]]
  s <- legacy_to_study(as.data.frame(get(e$tab, env)[[5]]), list(id = k, feature_type = "otu", subject = "personid"))$study
  a <- fit_study(s, e$f, rank = "phylum", longitudinal = e$long, n_cycles = 20)
  b <- fit_study(s, e$f, rank = "phylum", longitudinal = e$long, n_cycles = 500)
  ia <- attr(a, "fit_info"); ib <- attr(b, "fit_info")
  a <- a[a$term == e$term & a$feature %in% phyla, ]; b <- b[b$term == e$term & b$feature %in% phyla, ]
  m <- merge(a[, c("feature", "estimate", "se", "p_value", "converged")], b[, c("feature", "estimate", "se", "p_value", "converged")],
             by = "feature", suffixes = c("_20", "_500"))
  m$example <- k; m$d_est <- m$estimate_500 - m$estimate_20; m$d_p <- m$p_value_500 - m$p_value_20
  m$sig_20 <- m$p_value_20 < 0.05; m$sig_500 <- m$p_value_500 < 0.05
  out[[k]] <- m
}
r <- do.call(rbind, out); r$feature <- sub("^k__bacteria\\.p__", "", r$feature)
cat("\nMain phyla: converged at 20 cycles (the 1.x default) vs 500 cycles\n")
print(r[, c("example", "feature", "converged_20", "converged_500", "estimate_20", "estimate_500", "p_value_20", "p_value_500")], row.names = FALSE, digits = 3)
cat(sprintf("\nRows reported in the paper from fits NOT converged at 20 cycles: %d of %d\n", sum(!r$converged_20), nrow(r)))
cat(sprintf("Of those, converged with 500 cycles: %d\n", sum(!r$converged_20 & r$converged_500)))
cat(sprintf("Max |change in estimate| over all 16 rows: %.4f; max |change in p|: %.4f\n", max(abs(r$d_est)), max(abs(r$d_p))))
cat(sprintf("Rows whose significance at 0.05 changes: %d\n", sum(r$sig_20 != r$sig_500)))
