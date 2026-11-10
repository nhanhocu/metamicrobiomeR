suppressMessages(library(metamicrobiomeR))
S <- Sys.getenv("S"); st <- readRDS(file.path(S, "crc", "crc_studies_genus.rds"))
q <- function(expr) withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"), message = function(m) invokeRestart("muffleMessage"))
## A. Zeller: the same people on 16S and shotgun
z16 <- st$ZellerG_2014_16S; zsg <- st$ZellerG_2014
sg <- readRDS(file.path(S, "crc", "shotgun_plain.rds"))$ZellerG_2014$samples
subj_sg <- setNames(sg$subject_id, rownames(sg))
shared <- intersect(rownames(mm_samples(z16)), subj_sg[rownames(mm_samples(zsg))])
z16s <- subset_samples(z16, rownames(mm_samples(z16)) %in% shared)
zsgs <- subset_samples(zsg, subj_sg[rownames(mm_samples(zsg))] %in% shared)
cat("Zeller shared people:", length(shared), " 16S:", ncol(mm_counts(z16s)), " shotgun:", ncol(mm_counts(zsgs)), "\n")
print(table(mm_samples(z16s)$group)); print(table(mm_samples(zsgs)$group))
zres <- list()
for (e in c("lm_clr", "lm_asin", "bezi_gamlss")) {
  a <- as.data.frame(q(fit_study(z16s, ~ group, engine = e))); b <- as.data.frame(q(fit_study(zsgs, ~ group, engine = e)))
  m <- merge(a[, c("feature", "estimate", "se", "p_value")], b[, c("feature", "estimate", "se", "p_value")], by = "feature", suffixes = c("_16S", "_shotgun"))
  m$engine <- e; zres[[e]] <- m
  cat(sprintf("%-12s shared genera %3d  Pearson r %.2f  Spearman r %.2f  same sign %.0f%%  both p<0.05 %d, 16S only %d, shotgun only %d\n", e, nrow(m),
      cor(m$estimate_16S, m$estimate_shotgun), cor(m$estimate_16S, m$estimate_shotgun, method = "spearman"),
      100 * mean(sign(m$estimate_16S) == sign(m$estimate_shotgun)), sum(m$p_value_16S < 0.05 & m$p_value_shotgun < 0.05),
      sum(m$p_value_16S < 0.05 & m$p_value_shotgun >= 0.05), sum(m$p_value_16S >= 0.05 & m$p_value_shotgun < 0.05)))
}
saveRDS(zres, file.path(S, "crc", "zeller_platform.rds"))
## B. One-stage analysis of merged tables, against the two-stage result
main <- st[names(st) != "ZellerG_2014_16S"]
mres <- list()
for (fs in c("intersection", "union")) {
  mg <- q(merge_studies(main, features = fs))
  cat("merged (", fs, "):", nrow(mm_counts(mg)) - 1, "genera,", ncol(mm_counts(mg)), "samples\n")
  for (e in c("lm_clr", "lm_asin")) {
    r <- as.data.frame(q(fit_study(mg, ~ group + study, engine = e)))
    r <- r[r$term == "groupCRC", ]; r$p_adjust <- p.adjust(r$p_value, "fdr"); r$features <- fs; r$engine <- e; mres[[paste(fs, e)]] <- r
    cat(sprintf("  %-8s genera fitted %3d  FDR<0.05 %3d (up %d, down %d)\n", e, nrow(r), sum(r$p_adjust < 0.05), sum(r$p_adjust < 0.05 & r$estimate > 0), sum(r$p_adjust < 0.05 & r$estimate < 0)))
  }
}
saveRDS(mres, file.path(S, "crc", "merged_onestage.rds"))
