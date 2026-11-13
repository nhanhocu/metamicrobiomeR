# Pilot benchmark of the per-study engines on simulated data with known truth.
# Usage: Rscript data-raw/benchmark-pilot.R <engine> <scenario> <reps>
#   engines:   lm_asin | lm_clr | bezi_gamlss | zibb_glmmtmb
#   scenarios: null | alt_balanced | alt_unbalanced
# Results are cached in dev/cache/bench/; summarise with data-raw/benchmark-summary.R.
suppressMessages(library(metamicrobiomeR))
args <- commandArgs(trailingOnly = TRUE)
engine <- args[1]; scenario <- args[2]; reps <- as.integer(args[3])
sc <- list(null = list(n_diff = 0, effect = 0, balanced = TRUE),
           alt_balanced = list(n_diff = 4, effect = 1, balanced = TRUE),
           alt_unbalanced = list(n_diff = 4, effect = 1.5, balanced = FALSE))[[scenario]]
set.seed(100)
res <- vector("list", reps); t_all <- 0
for (r in seq_len(reps)) {
  sim <- simulate_study(n_samples = 60, n_features = 12, n_diff = sc$n_diff, effect = sc$effect, balanced = sc$balanced,
                        depth = 20000, depth_sd = 0.8, overdispersion = 20, seed = 1000 * match(scenario, c("null", "alt_balanced", "alt_unbalanced")) + r)
  t0 <- Sys.time()
  e <- suppressWarnings(fit_study(sim$study, ~ group, engine = engine, longitudinal = FALSE))
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  e <- e[e$term == "groupb" & e$part == "abundance", ]
  fi <- attr(e, "fit_info")
  d <- data.frame(rep = r, feature = e$feature, p = e$p_value, estimate = e$estimate, converged = e$converged, seconds = sec,
                  diff = sim$truth$diff[match(e$feature, sim$truth$feature)], n_features_fitted = nrow(e))
  res[[r]] <- d
}
out <- do.call(rbind, res); out$engine <- engine; out$scenario <- scenario
dir.create("dev/cache/bench", showWarnings = FALSE, recursive = TRUE)
saveRDS(out, file.path("dev/cache/bench", paste0(engine, "-", scenario, ".rds")))
cat(sprintf("%s / %s: %d reps, %d fits, %.1f s per replicate\n", engine, scenario, reps, nrow(out), mean(tapply(out$seconds, out$rep, mean))))
