# Q4 of the benchmark protocol, adapted from data-raw/benchmark-crossplatform.R (Payton Yau). Usage: Rscript benchmark-q4-crossplatform.R <replicates> <outdir>
# Restartable: replicates already in the output file are skipped.
# Cross-platform benchmark: two-stage meta-analysis versus a combined analysis (merged tables), repeated simulations.
# Each replicate: 3 studies on different platforms (genus lists, naming, depth, platform distortion), person-to-person variation,
# 5 truly different genera present in all studies, plus one genus present in two studies only.
suppressMessages(library(metamicrobiomeR)); options(metamicrobiomeR.quiet = TRUE)
a <- commandArgs(TRUE); R <- as.integer(a[1]); outdir <- a[2]; dir.create(outdir, showWarnings = FALSE, recursive = TRUE); out <- file.path(outdir, "Q4-crossplatform.csv")
done_seeds <- if (file.exists(out)) unique(read.csv(out)$seed) else integer(0)
G <- sprintf("Genus%02d", 1:60)
eff <- setNames(rep(0, 60), G); eff[G[c(20, 25, 40)]] <- 1.5; eff[G[c(30, 35)]] <- -1.3; eff[G[3]] <- 1.2
shared_true <- G[c(20, 25, 30, 35, 40)]
one_rep <- function(seed) {
  set.seed(seed); base <- rnorm(60, -4.5, 1.4)
  mk <- function(id, genera, depth, platform, decorate, n = 50, bias_sd = 0.4, person_sd = 1.0) {
    bias <- rnorm(length(genera), 0, bias_sd); grp <- factor(rep(c("a", "b"), each = n / 2)); cn <- matrix(0L, length(genera), n)
    for (j in seq_len(n)) { lp <- base[match(genera, G)] + bias + rnorm(length(genera), 0, person_sd) + (grp[j] == "b") * eff[match(genera, G)]
      p <- exp(lp) / sum(exp(lp)); cn[, j] <- rmultinom(1, round(depth * exp(rnorm(1, 0, 0.35))), p) }
    rownames(cn) <- paste0(id, "_f", seq_along(genera)); colnames(cn) <- paste0("S", seq_len(n))
    as_study(counts = cn, taxonomy = data.frame(genus = decorate(genera), row.names = rownames(cn)), samples = data.frame(group = grp, row.names = colnames(cn)),
             study = list(id = id, feature_type = "taxon", platform = platform))
  }
  studies <- list(mk("A", G[1:45], 3e4, "MiSeq", function(x) { x[5] <- paste0("[", x[5], "]"); x }), mk("B", G[1:50], 2e4, "IonTorrent", identity),
                  mk("C", G[16:60], 5e5, "shotgun", function(x) paste0(x, ifelse(seq_along(x) %% 3 == 0, "_A", ""))))
  rows <- list()
  score <- function(feat, est, padj) {
    nulls <- !feat %in% c(shared_true, G[3], G[55]); tr <- feat %in% shared_true
    if (sum(tr) == 0) return(c(auc = NA, dir = NA, power = NA, false_pos = NA))
    a <- mean(outer(abs(est[tr]), abs(est[nulls]), ">")) + 0.5 * mean(outer(abs(est[tr]), abs(est[nulls]), "=="))
    truth_sign <- sign(eff[match(feat[tr], G)])
    c(auc = a, dir = mean(sign(est[tr]) == truth_sign), power = sum(padj[tr] <= 0.05, na.rm = TRUE) / 5, false_pos = mean(padj[nulls] <= 0.05, na.rm = TRUE))
  }
  for (engine in c("zibb_glmmtmb", "lm_clr")) {
    t0 <- Sys.time()
    res <- suppressWarnings(meta_analyse(studies, ~ group, term = "groupb", rank = "genus", engine = engine, fit_args = list(longitudinal = FALSE), quality = "off", check = FALSE))
    p <- as.data.frame(res); sa <- score(p$feature, p$estimate, p$p_adjust)
    pz <- as.data.frame(pool_effects(res$effects, term = "groupb", method = "DL", test = "z")); sz <- score(pz$feature, pz$estimate, pz$p_adjust)
    rows[[length(rows) + 1]] <- data.frame(seed = seed, engine = engine, design = "two-stage (default test)", n_analysed = nrow(p), t(sa), seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
    rows[[length(rows) + 1]] <- data.frame(seed = seed, engine = engine, design = "two-stage (z test)", n_analysed = nrow(pz), t(sz), seconds = NA)
    t0 <- Sys.time()
    m <- suppressMessages(merge_studies(studies, rank = "genus"))
    e1 <- suppressWarnings(fit_study(m, ~ group + study, engine = engine, longitudinal = FALSE)); e1 <- e1[e1$term == "groupb", ]
    s1 <- score(e1$feature, e1$estimate, p.adjust(e1$p_value, "fdr"))
    rows[[length(rows) + 1]] <- data.frame(seed = seed, engine = engine, design = "combined (merged tables)", n_analysed = nrow(e1), t(s1), seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }
  do.call(rbind, rows)
}
for (i in seq_len(R)) { if ((1000 + i) %in% done_seeds) next; r <- try(one_rep(1000 + i), silent = TRUE); if (!inherits(r, "try-error")) write.table(r, out, sep = ",", row.names = FALSE, col.names = !file.exists(out), append = file.exists(out)) }
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo-q4.txt")); cat("all done\n")
