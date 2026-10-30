suppressMessages(library(metamicrobiomeR)); options(metamicrobiomeR.quiet = TRUE)
base <- simulate_studies(n_studies = 5, n_samples = 80, n_features = 40, n_diff = 6, effect = 1.2, tau = 0.0, depth_sd = 0.3, seed = 21)
truth <- base$truth$feature[base$truth$diff]
distort <- function(mode) lapply(seq_along(base$studies), function(i) {
  s <- base$studies[[i]]; set.seed(100 + i); f <- exp(rnorm(nrow(mm_counts(s)), 0, 1.0)); cn <- mm_counts(s); g <- mm_samples(s)$group
  if (mode == "study-wide") cn <- sweep(cn, 1, f, "*")                                    # same distortion for every sample of the study
  if (mode == "differs by group" && i >= 3) cn[, g == "b"] <- sweep(cn[, g == "b"], 1, f, "*")   # distortion applies to one group only, in 3 of 5 studies
  as_study(counts = round(cn), samples = mm_samples(s), study = list(id = s$study$id, feature_type = "asv")) })
cat(sprintf("%-24s %9s %12s %14s %10s\n", "batch effect", "median I2", "Q p < 0.05", "found (FDR)", "true found"))
for (mode in c("none", "study-wide", "differs by group")) {
  r <- meta_analyse(distort(mode), ~ group, term = "groupb", engine = "lm_asin", test = "z", quality = "off", check = FALSE, fit_args = list(longitudinal = FALSE, min_prevalence = 0))
  p <- as.data.frame(r); hit <- p$feature[p$p_adjust <= 0.05]
  cat(sprintf("%-24s %8.0f%% %11.0f%% %14d %7d of 6\n", mode, median(p$I2), 100 * mean(p$Q_p < 0.05), length(hit), sum(hit %in% truth)))
}
