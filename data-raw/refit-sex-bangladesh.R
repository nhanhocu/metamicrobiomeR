# Hybrid end-to-end check of the sex-difference meta-analysis.
# The sample data (gender, feeding, subject) of Haiti, USA CA/FL and USA NC are not in the
# repository, so only Bangladesh can be refitted from raw data. Its estimates are combined
# with the stored per-study results of the other three studies and pooled.
source("data-raw/legacy-helpers.R")
env <- new.env(); for (f in list.files("data-raw/legacy-1.1/results", full.names = TRUE)) load(f, env)
data(tabsex4); tab <- as.data.frame(tabsex4)

s <- legacy_to_study(as.data.frame(get("taxtab6.rm", env)[[5]]), list(id = "Bangladesh", feature_type = "otu", subject = "personid"))$study
dir.create("dev/validation", showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
fit <- fit_study(s, ~ gender + bf + age.sample, longitudinal = TRUE, n_cycles = 20, exclude_unassigned = FALSE)
cat(sprintf("Bangladesh refit: %d taxa, %.0f s\n", length(unique(fit$feature)), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
saveRDS(fit, "dev/validation/bangladesh-sex-fit.rds")

# (1) against the stored per-study result
st <- get("taxacom6.zi.rm.sex.adjustbfage", env)
cat("same taxa as stored:", setequal(unique(fit$feature), st$id), "\n")
for (tm in c("genderMale", "bfNon_exclusiveBF", "bfNo_BF", "age.sample")) {
  x <- fit[fit$term == tm, ]; x <- x[match(st$id, x$feature), ]
  cat(sprintf("  %-18s max |diff| estimate %.1e, se %.1e, p %.1e\n", tm, max(abs(x$estimate - st[[paste0("Estimate.", tm)]]), na.rm = TRUE),
              max(abs(x$se - st[[paste0("Std. Error.", tm)]]), na.rm = TRUE), max(abs(x$p_value - st[[paste0("Pr(>|t|).", tm)]]), na.rm = TRUE)))
}

# (2) pool: refitted Bangladesh + stored results of the other three studies
p <- parse_legacy_names(tab$id); rank_of <- unname(codes)[p$level]
est <- tab[["Estimate.genderMale"]]; se <- tab[["Std. Error.genderMale"]]; ok <- !is.na(est) & !is.na(se) & tab$pop != "Bangladesh"
others <- new_effect_table(study = tab$pop[ok], feature = tab$id[ok], term = "genderMale", estimate = est[ok], se = se[ok], rank = rank_of[ok],
                           part = "abundance", engine = "bezi_gamlss", effect_scale = "logit_mean")
mine <- fit[fit$term == "genderMale", ]
effects <- bind_effects(mine, others)
pooled <- pool_effects(effects, method = "DL", test = "z", adjust_by = c("part", "term"))
old <- get("metab.sex", env)$random$genderMale
m <- merge(pooled, old, by.x = "feature", by.y = "id", suffixes = c("", "_old"))
cat(sprintf("\nPooled (refitted Bangladesh + 3 stored studies): %d taxa vs %d stored; same taxa: %s\n", nrow(pooled), nrow(old), setequal(pooled$feature, old$id)))
cat(sprintf("  max |diff| estimate %.1e, se %.1e, p %.1e, p.adjust %.1e\n", max(abs(m$estimate - m$estimate_old)), max(abs(m$se - m$se_old)),
            max(abs(m$p_value - m$p)), max(abs(m$p_adjust - m$p.adjust))))
