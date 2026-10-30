# Mock-data example of a full meta-analysis with metamicrobiomeR 2.0
# Four simulated studies; 5 of 30 features are truly different between groups (effect size 1.5 on the logit scale, between-study SD 0.3).
suppressMessages(library(metamicrobiomeR))
sim <- simulate_studies(n_studies = 4, n_samples = 60, n_features = 30, n_diff = 5, effect = 1.5, tau = 0.3, depth_sd = 0.4, seed = 2026)
studies <- sim$studies
cat("TRUTH: differential features are", paste(sim$truth$feature[sim$truth$diff], collapse = ", "), "\n\n")

cat("---- 1. Are the studies safe to pool? ----\n")
print(check_poolable(studies))

cat("\n---- 2. Differential abundance: one call ----\n")
res <- meta_analyse(studies, ~ group, term = "groupb", engine = "lm_asin", fit_args = list(longitudinal = FALSE))
print(res)
tab <- as.data.frame(res); hit <- tab$feature[tab$p_adjust <= 0.05]
truth <- sim$truth$feature[sim$truth$diff]
cat(sprintf("\nFound %d at adjusted p <= 0.05: %s\nTrue positives: %d of 5 | false positives: %d\n", length(hit), paste(hit, collapse = ", "), sum(hit %in% truth), sum(!hit %in% truth)))

cat("\n---- 3. Does the answer depend on the test? ----\n")
print(pool_sensitivity(res$effects, term = "groupb", significance = "p_adjust"))

cat("\n---- 4. Alpha diversity ----\n")
alpha <- bind_effects(lapply(studies, fit_alpha, ~ group, longitudinal = FALSE))
print(pool_effects(alpha, term = "groupb", min_prop = 0))

cat("\n---- 5. Beta diversity ----\n")
if (requireNamespace("vegan", quietly = TRUE)) {
  beta <- lapply(studies, fit_beta, ~ group, distance = "bray", permutations = 499, seed = 1)
  print(do.call(pool_beta, beta))
}
cat("\n---- 6. Moderator: does the effect differ by study? ----\n")
sd_ <- data.frame(study = res$studies$study, depth = c(1, 1, 2, 2))
print(meta_regress(res$effects, sd_, ~ depth, term = "groupb", test = "z")[1:3, c("feature", "moderator", "estimate", "p_value", "QM_p")])

dir.create("out", showWarnings = FALSE)
ggplot2::ggsave("out/forest.png", plot(res, top = 3), width = 7, height = 7, dpi = 110, bg = "white")
ggplot2::ggsave("out/heatmap.png", plot(res, "heatmap", top = 12), width = 7, height = 6, dpi = 110, bg = "white")
ggplot2::ggsave("out/volcano.png", plot(res, "volcano"), width = 6.5, height = 5, dpi = 110, bg = "white")
cat("\nsaved plots\n")
