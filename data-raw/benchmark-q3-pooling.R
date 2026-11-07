# Benchmark Q3 (protocol, sections 2 to 8): calibration and power of the pooled tests with 3, 5 and 10 studies.
#
# Usage: Rscript data-raw/benchmark-q3-pooling.R <replicates> <output-dir> [cores]
#
# Studies come from simulate_studies() (generator A, counts) with between-study heterogeneity tau on the effect. Each study
# is fitted with lm_asin, the fastest engine, so that the pooling step is what varies. With effect = 0 the mean of the
# study effects is 0, so a pooled p-value below 0.05 is a false positive whatever tau is.

suppressPackageStartupMessages({ library(metamicrobiomeR); library(parallel) })
options(metamicrobiomeR.quiet = TRUE)
a <- commandArgs(TRUE)
nrep <- as.integer(a[1]); out_dir <- a[2]; cores <- if (length(a) >= 3) as.integer(a[3]) else 2L
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

grid <- expand.grid(k = c(3, 5, 10), tau = c(0, 0.3), effect = c(0, 1), stringsAsFactors = FALSE)
grid$scenario <- sprintf("Q3-k%02d-tau%02d-eff%02d", grid$k, round(10 * grid$tau), round(10 * grid$effect))
methods <- list(REML_adhoc = list(method = "REML", test = "adhoc"), REML_z = list(method = "REML", test = "z"),
                DL_z = list(method = "DL", test = "z"))

one <- function(r, sc) {
  sim <- simulate_studies(n_studies = sc$k, tau = sc$tau, n_samples = 100, n_features = 20, n_diff = if (sc$effect == 0) 0 else 2,
                          effect = sc$effect, seed = 30000 + r)
  fits <- lapply(sim$studies, function(st) suppressWarnings(suppressMessages(
    fit_study(st, ~ group, engine = "lm_asin", longitudinal = FALSE, min_prevalence = 0.1))))
  eff <- suppressMessages(do.call(bind_effects, fits))
  do.call(rbind, lapply(names(methods), function(nm) {
    pl <- suppressWarnings(suppressMessages(tryCatch(
      pool_effects(eff, term = "groupb", part = "abundance", method = methods[[nm]]$method, test = methods[[nm]]$test),
      error = function(e) NULL)))
    if (is.null(pl) || !nrow(pl)) return(NULL)
    m <- match(pl$feature, sim$truth$feature)
    data.frame(rep = r, method = nm, feature = pl$feature, k_used = pl$k, estimate = pl$estimate, se = pl$se,
               ci_lo = pl$ci_lo, ci_hi = pl$ci_hi, p = pl$p_value, diff = sim$truth$diff[m], beta = sim$truth$beta[m])
  }))
}

only <- Sys.getenv("BENCH_ONLY")
if (nzchar(only)) grid <- grid[grid$scenario %in% strsplit(only, ",")[[1]], , drop = FALSE]
for (i in seq_len(nrow(grid))) {
  sc <- grid[i, ]
  f <- file.path(out_dir, paste0(sc$scenario, ".csv"))
  if (file.exists(f)) next
  res <- do.call(rbind, mclapply(seq_len(nrep), one, sc = sc, mc.cores = cores))
  res <- cbind(scenario = sc$scenario, res, sc[c("k", "tau", "effect")], row.names = NULL)
  utils::write.csv(res, f, row.names = FALSE)
  cat(sc$scenario, "done", format(Sys.time()), "\n")
}
writeLines(utils::capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo-q3.txt"))
cat("all done\n")
