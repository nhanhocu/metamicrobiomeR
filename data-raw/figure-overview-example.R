# Example overview plots for the 2019 sex-difference meta-analysis (genus level, stored per-study results)
source("data-raw/legacy-helpers.R")
load("data-raw/legacy-cran-1.2-data/tabsex4.rda"); tab <- as.data.frame(tabsex4); p <- parse_legacy_names(tab$id)
est <- tab[["Estimate.genderMale"]]; se <- tab[["Std. Error.genderMale"]]; ok <- !is.na(est) & !is.na(se) & unname(codes)[p$level] == "genus"
effects <- new_effect_table(study = tab$pop[ok], feature = tab$id[ok], term = "genderMale", estimate = est[ok], se = se[ok],
                            rank = "genus", part = "abundance", p_value = tab[["Pr(>|t|).genderMale"]][ok], engine = "bezi_gamlss", effect_scale = "logit_mean")
short <- function(x) sub(".*g__", "", x); lab <- stats::setNames(short(unique(effects$feature)), unique(effects$feature))
pooled <- pool_effects(effects)
dir.create("dev/figures", showWarnings = FALSE)
h <- plot_heatmap(effects, pooled, top = 15, feature_labels = lab)
ggplot2::ggsave("dev/figures/heatmap-sex-example.png", h, width = 7.5, height = 6, dpi = 130, bg = "white")
v <- plot_volcano(pooled, feature_labels = lab)
ggplot2::ggsave("dev/figures/volcano-sex-example.png", v, width = 6.5, height = 5, dpi = 130, bg = "white")
cat("saved\n")
