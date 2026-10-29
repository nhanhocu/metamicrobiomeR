# Reproduce the 2019 sex-difference meta-analysis (Fig. 5) from the stored per-study results.
# Part A (this script): pool the stored per-study estimates with pool_effects() and compare
# with the stored pooled result (`metab.sex`, produced by meta.taxa() in 1.x).
source("data-raw/legacy-helpers.R")
env <- new.env(); for (f in list.files("data-raw/legacy-1.1/results", full.names = TRUE)) load(f, env)
load("data-raw/legacy-cran-1.2-data/tabsex4.rda")      # the four per-study results, stacked (bundled with CRAN 1.2)
tab <- as.data.frame(tabsex4)
terms <- sub("^Estimate\\.", "", grep("^Estimate\\.", colnames(tab), value = TRUE))
p <- parse_legacy_names(tab$id); rank_of <- unname(codes)[p$level]

to_effects <- function(term) {
  est <- tab[[paste0("Estimate.", term)]]; se <- tab[[paste0("Std. Error.", term)]]
  ok <- !is.na(est) & !is.na(se)
  new_effect_table(study = tab$pop[ok], feature = tab$id[ok], term = term, estimate = est[ok], se = se[ok],
                   rank = rank_of[ok], part = "abundance", p_value = tab[[paste0("Pr(>|t|).", term)]][ok],
                   engine = "bezi_gamlss", effect_scale = "logit_mean")
}
effects <- bind_effects(lapply(terms, to_effects))
cat("studies:", paste(unique(effects$study), collapse = ", "), "| rows:", nrow(effects), "\n")

# 1.x settings: DerSimonian-Laird, normal-theory CI, taxa in >= 50% of studies, FDR over all taxa of all ranks
pooled <- pool_effects(effects, method = "DL", test = "z", min_prop = 0.5, adjust_by = c("part", "term"))
cat(sprintf("\n%-20s %7s %7s | max |difference| vs stored 1.x result (estimate, se, ll, ul, p, p.adjust)\n", "term", "new n", "stored"))
for (tm in terms) {
  s <- get("metab.sex", env)$random[[tm]]; n <- pooled[pooled$term == tm, ]
  m <- merge(n, s, by.x = "feature", by.y = "id", suffixes = c("", "_old"))
  d <- c(estimate = max(abs(m$estimate - m$estimate_old)), se = max(abs(m$se - m$se_old)), ll = max(abs(m$ci_lo - m$ll)),
         ul = max(abs(m$ci_hi - m$ul)), p = max(abs(m$p_value - m$p)), p.adjust = max(abs(m$p_adjust - m$p.adjust)))
  cat(sprintf("%-20s %7d %7d | %s | same taxa: %s\n", tm, nrow(n), nrow(s), paste(format(d, digits = 2, scientific = TRUE), collapse = "  "),
              setequal(n$feature, s$id)))
}

# --- Why do the other terms differ in the set of pooled taxa? -------------------------------
cat("\n== taxa pooled by 1.x from fewer than 2 studies with an estimate for the term ==\n")
for (tm in terms[-1]) {
  s <- get("metab.sex", env)$random[[tm]]
  e <- effects[effects$term == tm, ]
  k <- table(e$feature)[s$id]; k[is.na(k)] <- 0
  cat(sprintf("%-20s stored pooled taxa: %3d | with 0 estimates: %3d | with exactly 1: %3d | with >= 2: %3d\n", tm, nrow(s), sum(k == 0), sum(k == 1), sum(k >= 2)))
}

# --- The paper's headline conclusions (sex difference, 4 studies) ---------------------------
g <- pooled[pooled$term == "genderMale", ]
cat("\n== genderMale, 1.x settings (DL, normal CI, FDR across all ranks) ==\n")
cat("phylum-level taxa with pooled p < 0.05:", sum(g$rank == "phylum" & g$p_value < 0.05), "of", sum(g$rank == "phylum"), "\n")
gg <- g[g$rank == "genus", ]
cat("genera with pooled p < 0.05:", sum(gg$p_value < 0.05), "of", nrow(gg), "->", paste(sub(".*\\.g__", "", gg$feature[gg$p_value < 0.05]), collapse = ", "), "\n")
cat("genera with FDR < 0.05:", paste(sub(".*\\.g__", "", gg$feature[gg$p_adjust < 0.05]), collapse = ", "),
    "| smallest FDR:", format(min(gg$p_adjust), digits = 2), "\n")

# --- Same data with the 2.0 defaults (REML, Hartung-Knapp, FDR within rank) ---------------------
n2 <- pool_effects(effects, term = "genderMale")
cat("\n== genderMale, 2.0 defaults (REML, ad hoc Hartung-Knapp, FDR within rank, unassigned lineages kept) ==\n")
cat("phylum-level taxa with p < 0.05:", sum(n2$rank == "phylum" & n2$p_value < 0.05), "\n")
ng <- n2[n2$rank == "genus", ]
cat("genera with p < 0.05:", sum(ng$p_value < 0.05), "->", paste(sub(".*\\.g__", "", ng$feature[ng$p_value < 0.05]), collapse = ", "), "\n")
cat("genera with FDR < 0.05:", paste(sub(".*\\.g__", "", ng$feature[ng$p_adjust < 0.05]), collapse = ", "), "\n")
cop <- ng[grepl("coprococcus", ng$feature), ]
cat(sprintf("Coprococcus: DL/z p = %s, FDR = %s | 2.0 p = %s, FDR = %s; I2 = %.0f%%\n",
            format(gg$p_value[grepl("coprococcus", gg$feature)], digits = 2), format(gg$p_adjust[grepl("coprococcus", gg$feature)], digits = 2),
            format(cop$p_value, digits = 2), format(cop$p_adjust, digits = 2), cop$I2))
