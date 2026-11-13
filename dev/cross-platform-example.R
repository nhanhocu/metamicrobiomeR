suppressMessages(library(metamicrobiomeR))
set.seed(42)
G <- sprintf("Genus%02d", 1:60); diff_shared <- G[c(20, 25, 30, 35, 40)]; diff_unshared <- G[c(3, 55)]   # 5 truly different genera in all studies, 2 where only some studies have them
effect <- setNames(rep(0, 60), G); effect[c(diff_shared, diff_unshared)] <- 1.2
base <- rnorm(60, -4, 1.5)
make_study <- function(id, genera, depth, platform, region, db, decorate, n = 40) {
  bias <- rnorm(length(genera), 0, 0.5)                                       # platform-specific distortion, different for every genus
  grp <- factor(rep(c("a", "b"), each = n / 2)); cn <- matrix(0L, length(genera), n)
  for (j in seq_len(n)) {
    lp <- base[match(genera, G)] + bias + (grp[j] == "b") * effect[match(genera, G)]
    p <- exp(lp) / sum(exp(lp)); cn[, j] <- rmultinom(1, round(depth * exp(rnorm(1, 0, 0.3))), p)
  }
  nm <- decorate(genera)
  rownames(cn) <- paste0(id, "_f", seq_along(genera)); colnames(cn) <- paste0("S", seq_len(n))
  as_study(counts = cn, taxonomy = data.frame(genus = nm, row.names = rownames(cn)), samples = data.frame(group = grp, row.names = colnames(cn)),
           study = list(id = id, feature_type = "taxon", platform = platform, region = region, database = db))
}
A <- make_study("A", G[1:45], 2e4, "Illumina MiSeq", "V4", "SILVA 138", function(x) { x[5] <- paste0("[", x[5], "]"); x })
B <- make_study("B", G[1:50], 3e4, "Ion Torrent", "V3-V4", "SILVA 138", function(x) x)
C <- make_study("C", G[16:60], 2e6, "NovaSeq shotgun", NA, "GTDB r207", function(x) paste0(x, ifelse(seq_along(x) %% 3 == 0, "_A", "")))
studies <- list(A, B, C)
truth_all <- c(diff_shared, diff_unshared)
cat("What differs between the studies:\n"); print(study_table(studies)[, c("study", "n_samples", "n_features", "platform", "region", "database")], row.names = FALSE)
cat("\n== Can they be pooled? ==\n"); print(check_poolable(studies, rank = "genus"))

auc <- function(est, feat) { pos <- est[feat %in% diff_shared]; neg <- est[!feat %in% truth_all]; mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "==")) }   # P(a true genus has a larger estimate than a null genus)
score <- function(label, feat, est, padj, note = "") {
  cat(sprintf("%-30s analysed %2d | effect AUC %.2f | median estimate: true %+.2f, null %+.2f | FDR 5%%: %2d found%s\n", label, length(feat), auc(est, feat),
              median(est[feat %in% diff_shared]), median(est[!feat %in% truth_all]), sum(padj <= 0.05, na.rm = TRUE), note)) }
cat("\n== Does the engine matter across platforms? (two-stage, 3 studies) ==\n")
for (eng in c("lm_asin", "zibb_glmmtmb")) {
  res <- suppressWarnings(suppressMessages(meta_analyse(studies, ~ group, term = "groupb", rank = "genus", engine = eng, fit_args = list(longitudinal = FALSE), quality = "off", check = FALSE, test = "z")))
  p2 <- as.data.frame(res)
  score(sprintf("two-stage, %s", eng), p2$feature, p2$estimate, p2$p_adjust, sprintf(" | median I2 %.0f%%", median(p2$I2)))
}
cat("\n== Combined analysis (merge the raw tables, study as a covariate) ==\n")
m <- suppressMessages(merge_studies(studies, rank = "genus"))
for (eng in c("lm_asin", "zibb_glmmtmb")) {
  e1 <- suppressWarnings(fit_study(m, ~ group + study, engine = eng, longitudinal = FALSE)); e1 <- e1[e1$term == "groupb", ]
  score(sprintf("combined, %s", eng), e1$feature, e1$estimate, p.adjust(e1$p_value, "fdr"))
}
cat("\nGenus03 (in the two 16S studies only): two-stage can pool it from 2 studies; the combined analysis drops it (not in every study).\n")
