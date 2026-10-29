# Do the new diagnostics catch the USA(NC) problem in the 2019 sex-difference meta-analysis?
source("data-raw/legacy-helpers.R")
load("data-raw/legacy-cran-1.2-data/tabsex4.rda"); tab <- as.data.frame(tabsex4); p <- parse_legacy_names(tab$id)
est <- tab[["Estimate.genderMale"]]; se <- tab[["Std. Error.genderMale"]]; ok <- !is.na(est) & !is.na(se)
effects <- new_effect_table(study = tab$pop[ok], feature = tab$id[ok], term = "genderMale", estimate = est[ok], se = se[ok],
                            rank = unname(codes)[p$level][ok], part = "abundance", p_value = tab[["Pr(>|t|).genderMale"]][ok],
                            engine = "bezi_gamlss", effect_scale = "logit_mean")
cop <- "k__bacteria.p__firmicutes.c__clostridia.o__clostridiales.f__lachnospiraceae.g__coprococcus"

cat("== check_effects_quality ==\n"); q <- check_effects_quality(effects); print(q)
print(q$studies[, c("study", "n", "se_min", "se_q", "reference_se_q", "flagged_se")], digits = 3, row.names = FALSE)
cat("\nflagged taxa in USA(NC):", sum(q$rows$flag_se), "| in the other studies:", sum(q$rows$flag_se[q$rows$study != "USA(NC)"]), "\n")
cat("Coprococcus flagged in:", paste(q$rows$study[q$rows$feature == cop & q$rows$flag_se], collapse = ", "), "\n")

cat("\n== pool_sensitivity (genus level) ==\n")
gen <- effects[effects$rank %in% "genus", ]
s <- pool_sensitivity(gen, significance = "p_adjust", adjust_by = c("part", "term", "rank"))
print(s); print(s$wide[grepl("coprococcus$", s$wide$feature), c("k", "p_DL_z", "p_REML_knha", "p_REML_adhoc", "n_settings_significant")], digits = 3)

cat("\n== filter_effects then pool ==\n")
f <- suppressMessages(filter_effects(effects)); cat("dropped:", nrow(attr(f, "dropped")), "of", nrow(effects), "estimates (all in", paste(unique(attr(f, "dropped")$study), collapse = ", "), ")\n")
for (lab in c("before", "after")) {
  e <- if (lab == "before") effects else f
  pp <- pool_effects(e, method = "DL", test = "z", adjust_by = c("part", "term"))
  g <- pp[pp$rank == "genus" & pp$p_adjust < 0.05, ]
  cat(sprintf("%-6s filtering: genera with FDR < 0.05 (1.x settings): %d%s | Coprococcus: %s\n", lab, nrow(g), if (nrow(g)) paste0(" (", paste(sub(".*g__", "", g$feature), collapse = ", "), ")") else "",
              if (cop %in% pp$feature) sprintf("k = %d, p = %.2g", pp$k[pp$feature == cop], pp$p_value[pp$feature == cop]) else "not pooled (fewer than 2 studies)"))
}

cat("\n== combine_pvalues (needs only p-values and signs; weights sqrt(n); n is unknown here so equal weights) ==\n")
cp <- combine_pvalues(effects, weights = "none", adjust_by = c("part", "term"))
print(cp[cp$feature == cop, c("k", "n_pos", "n_neg", "direction", "sign_test_p", "statistic", "p_value", "p_adjust")], digits = 3)
cp2 <- combine_pvalues(f, weights = "none", adjust_by = c("part", "term"))
cat("after filtering: Coprococcus", if (cop %in% cp2$feature) sprintf("k = %d, combined p = %.3g", cp2$k[cp2$feature == cop], cp2$p_value[cp2$feature == cop]) else "not combined", "\n")
