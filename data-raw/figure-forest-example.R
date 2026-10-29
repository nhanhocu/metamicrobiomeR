# Example forest plot: sex difference in infant gut microbiome, four studies (stored per-study results)
source("data-raw/legacy-helpers.R")
data(tabsex4); tab <- as.data.frame(tabsex4); p <- parse_legacy_names(tab$id)
est <- tab[["Estimate.genderMale"]]; se <- tab[["Std. Error.genderMale"]]; ok <- !is.na(est) & !is.na(se)
effects <- new_effect_table(study = tab$pop[ok], feature = tab$id[ok], term = "genderMale", estimate = est[ok], se = se[ok],
                            rank = unname(codes)[p$level][ok], part = "abundance", engine = "bezi_gamlss", effect_scale = "logit_mean")
pooled <- pool_effects(effects)                                   # 2.0 defaults
feat <- c("k__bacteria.p__firmicutes.c__clostridia.o__clostridiales.f__lachnospiraceae.g__coprococcus",
          "k__bacteria.p__firmicutes.c__clostridia.o__clostridiales.f__veillonellaceae.g__megamonas",
          "k__bacteria.p__bacteroidetes")
stopifnot(all(feat %in% effects$feature))
labels <- c("Coprococcus (genus)", "Megamonas (genus)", "Bacteroidetes (phylum)"); names(labels) <- feat
g <- plot_forest(effects, pooled, feature = feat, term = "genderMale", feature_labels = labels)
ggplot2::ggsave("dev/figures/forest-sex-example.png", g, width = 7.5, height = 8, dpi = 150, bg = "white")
cat("saved; features:", length(feat), "\n")
