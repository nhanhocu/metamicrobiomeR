suppressMessages(library(metamicrobiomeR))
S <- Sys.getenv("S"); od <- file.path(S, "crc", "pooled"); rd <- function(t) readRDS(file.path(od, paste0(t, ".rds")))
pc <- c("Fusobacterium", "Parvimonas", "Peptostreptococcus", "Porphyromonas", "Gemella")
cat("== Positive controls, all 15 studies\n")
for (t in c("clr_all_REML_adhoc", "clr_all_REML_z", "asin_all_REML_adhoc", "bezi_all_REML_adhoc", "clr_adj_all_REML_adhoc")) {
  p <- rd(t); p <- p[match(pc, p$feature), ]
  cat(t, "\n"); print(data.frame(genus = pc, k = p$k, est = round(p$estimate, 2), lo = round(p$ci_lo, 2), hi = round(p$ci_hi, 2), p_adj = signif(p$p_adjust, 2), I2 = round(p$I2), maxw = round(p$max_weight, 2)), row.names = FALSE)
}
cat("\n== Significant genera, lm_clr, all studies, REML HK (default)\n")
p <- rd("clr_all_REML_adhoc"); z <- rd("clr_all_REML_z"); s <- p[p$p_adjust < 0.05, ]; s <- s[order(-s$estimate), ]
s$p_adj_z <- z$p_adjust[match(s$feature, z$feature)]
print(data.frame(genus = s$feature, k = s$k, est = round(s$estimate, 2), p_adj = signif(s$p_adjust, 2), I2 = round(s$I2), pi_lo = round(s$pi_lo, 2), pi_hi = round(s$pi_hi, 2), maxw = round(s$max_weight, 2), top = s$top_study), row.names = FALSE)
cat("\nOnly with z:", paste(setdiff(z$feature[z$p_adjust < 0.05], s$feature), collapse = ", "), "\n")
sh <- rd("clr_shotgun_REML_adhoc"); a16 <- rd("clr_16S_REML_z")
m <- merge(sh[, c("feature", "estimate")], a16[, c("feature", "estimate")], by = "feature", suffixes = c("_sg", "_16"))
cat(sprintf("\nShotgun-only vs 16S-only pooled lm_clr estimates: %d genera, Pearson r %.2f, same sign %.0f%%\n", nrow(m), cor(m$estimate_sg, m$estimate_16), 100 * mean(sign(m$estimate_sg) == sign(m$estimate_16))))
sig <- s$feature; mm <- m[m$feature %in% sig, ]
cat(sprintf("  among the %d HK-significant genera present in both: same sign %d\n", nrow(mm), sum(sign(mm$estimate_sg) == sign(mm$estimate_16))))
