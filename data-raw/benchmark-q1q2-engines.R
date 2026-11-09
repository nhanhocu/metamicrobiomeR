# Benchmark Q1 and Q2 (protocol, sections 2 to 8): calibration and power of the per-study engines.
#
# Usage: Rscript data-raw/benchmark-q1q2-engines.R <set> <replicates> <output-dir> [cores]
#   set = "cross" (generator A, counts, cross-sectional, all five engines)
#         "long"  (generator B, zero-inflated beta with a subject random intercept, one feature per data set)
#         "longc" (generator C, not the BEZI model: logistic-normal binomial counts with a subject random intercept,
#                  zeros from sampling only, one feature per data set; added after Payton's review of the protocol)
#
# Truth (protocol section 6): a feature is differential when its generating effect is not zero. Non-differential features
# in a data set with differential ones can still shift through compositional effects, and count as false positives here.

suppressPackageStartupMessages({ library(metamicrobiomeR); library(parallel) })
options(metamicrobiomeR.quiet = TRUE)
a <- commandArgs(TRUE)
set <- a[1]; nrep <- as.integer(a[2]); out_dir <- a[3]; cores <- if (length(a) >= 4) as.integer(a[4]) else 2L
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

quiet_fit <- function(...) suppressWarnings(suppressMessages(tryCatch(fit_study(...), error = function(e) NULL)))

# ---- Generator A: counts, cross-sectional, five engines ------------------------------------------------------------
cross_grid <- rbind(
  data.frame(scenario = "A-null-n30",    n = 30,  zi = 0,   effect = 0),
  data.frame(scenario = "A-null-n100",   n = 100, zi = 0,   effect = 0),
  data.frame(scenario = "A-null-n300",   n = 300, zi = 0,   effect = 0),
  data.frame(scenario = "A-null-zi50",   n = 100, zi = 0.5, effect = 0),
  data.frame(scenario = "A-eff05-n100",  n = 100, zi = 0,   effect = 0.5),
  data.frame(scenario = "A-eff10-n100",  n = 100, zi = 0,   effect = 1.0),
  data.frame(scenario = "A-eff10-zi50",  n = 100, zi = 0.5, effect = 1.0)
)
cross_engines <- c("bezi_gamlss", "bezi_glmmtmb", "zibb_glmmtmb", "lm_asin", "lm_clr")

cross_one <- function(r, sc) {
  sim <- simulate_study(n_samples = sc$n, n_features = 20, n_diff = if (sc$effect == 0) 0 else 2, effect = sc$effect,
                        zero_inflation = sc$zi, seed = 10000 + r)
  st <- sim$study
  share0 <- mean(mm_counts(st) == 0)
  do.call(rbind, lapply(cross_engines, function(eng) {
    t0 <- proc.time()[["elapsed"]]
    e <- quiet_fit(st, ~ group, engine = eng, longitudinal = FALSE, min_prevalence = 0.1)
    secs <- proc.time()[["elapsed"]] - t0
    if (is.null(e) || !nrow(e)) return(data.frame(rep = r, engine = eng, feature = NA, estimate = NA, se = NA, p = NA, converged = NA,
                                                   diff = NA, beta = NA, share0 = share0, secs = secs))
    e <- e[e$part == "abundance" & e$term == "groupb", ]
    m <- match(e$feature, sim$truth$feature)
    data.frame(rep = r, engine = eng, feature = e$feature, estimate = e$estimate, se = e$se, p = e$p_value, converged = e$converged,
               diff = sim$truth$diff[m], beta = sim$truth$beta[m], share0 = share0, secs = secs)
  }))
}

# ---- Generator B: zero-inflated beta with a subject random intercept, one feature -------------------------------------
long_grid <- rbind(
  expand.grid(effect = 0, prev = c(1, 0.6, 0.15), term = c("between", "within"), stringsAsFactors = FALSE),
  expand.grid(effect = c(0.5, 1), prev = 0.6, term = c("between", "within"), stringsAsFactors = FALSE)
)
long_grid$scenario <- sprintf("B-%s-eff%02d-prev%03d", long_grid$term, round(10 * long_grid$effect), round(100 * long_grid$prev))
long_engines <- list(bezi_gamlss = list(engine = "bezi_gamlss", longitudinal = TRUE),
                     bezi_gamlss_norandom = list(engine = "bezi_gamlss", longitudinal = FALSE),
                     bezi_glmmtmb = list(engine = "bezi_glmmtmb", longitudinal = TRUE),
                     lm_asin = list(engine = "lm_asin", longitudinal = TRUE))

long_one <- function(r, sc) {
  set.seed(20000 + r)
  n <- 40; m <- 4; id <- rep(seq_len(n), each = m)
  g <- if (sc$term == "between") rep(rep(c("a", "b"), length.out = n), each = m) else rep(c("a", "b"), length.out = n * m)
  age <- stats::runif(n * m, 0, 6); u <- stats::rnorm(n, sd = 0.5)[id]
  mu <- stats::plogis(-3 + u + 0.1 * age + sc$effect * (g == "b"))
  y <- stats::rbeta(n * m, mu * 20, (1 - mu) * 20); y[stats::runif(n * m) > sc$prev] <- 0
  sam <- data.frame(subj = paste0("s", id), g = g, age = age, row.names = paste0("S", seq_len(n * m)))
  rel <- rbind(F1 = y, other = 1 - y); colnames(rel) <- rownames(sam)
  st <- as_study(relabund = rel, samples = sam, study = list(id = "sim", feature_type = "otu", subject = "subj"))
  do.call(rbind, lapply(names(long_engines), function(nm) {
    cfg <- long_engines[[nm]]
    t0 <- proc.time()[["elapsed"]]
    e <- quiet_fit(st, ~ g + age, engine = cfg$engine, longitudinal = cfg$longitudinal, min_prevalence = 0, min_mean = 0,
                   exclude_unassigned = FALSE)
    secs <- proc.time()[["elapsed"]] - t0
    row <- if (is.null(e)) NULL else e[e$feature == "F1" & e$part == "abundance" & e$term == "gb", ]
    if (is.null(row) || !nrow(row)) return(data.frame(rep = r, engine = nm, estimate = NA, se = NA, p = NA, converged = NA, secs = secs))
    data.frame(rep = r, engine = nm, estimate = row$estimate, se = row$se, p = row$p_value, converged = row$converged, secs = secs)
  }))
}

# ---- Generator C: logistic-normal counts with a subject random intercept, not the BEZI model ---------------------------
# Feature F1 against the rest of the community. Its log odds are b0 + subject effect (SD 0.5) + 0.1 * age + effect * (g == b)
# + sample-level noise (SD 1, logistic-normal rather than beta), and its count is binomial on a library size of about 5,000
# reads (log-normal, SD 0.5). Zeros come from sampling only. b0 = -3, -8.9 and -11.2 give about 100%, 60% and 15% non-zero.
longc_grid <- rbind(
  expand.grid(effect = 0, prev = c(1, 0.6, 0.15), term = c("between", "within"), stringsAsFactors = FALSE),
  expand.grid(effect = c(0.5, 1), prev = 0.6, term = c("between", "within"), stringsAsFactors = FALSE)
)
longc_grid$scenario <- sprintf("C-%s-eff%02d-prev%03d", longc_grid$term, round(10 * longc_grid$effect), round(100 * longc_grid$prev))
longc_engines <- c(long_engines, list(zibb_glmmtmb = list(engine = "zibb_glmmtmb", longitudinal = TRUE),
                                      lm_clr = list(engine = "lm_clr", longitudinal = TRUE)))
longc_b0 <- c("1" = -3, "0.6" = -8.9, "0.15" = -11.2)

longc_one <- function(r, sc) {
  set.seed(40000 + r)
  n <- 40; m <- 4; id <- rep(seq_len(n), each = m)
  g <- if (sc$term == "between") rep(rep(c("a", "b"), length.out = n), each = m) else rep(c("a", "b"), length.out = n * m)
  age <- stats::runif(n * m, 0, 6); u <- stats::rnorm(n, sd = 0.5)[id]
  depth <- round(exp(stats::rnorm(n * m, log(5000), 0.5)))
  p <- stats::plogis(longc_b0[[as.character(sc$prev)]] + u + 0.1 * age + sc$effect * (g == "b") + stats::rnorm(n * m, sd = 1))
  y <- stats::rbinom(n * m, depth, p)
  sam <- data.frame(subj = paste0("s", id), g = g, age = age, row.names = paste0("S", seq_len(n * m)))
  cn <- rbind(F1 = y, other = depth - y); colnames(cn) <- rownames(sam)
  st <- as_study(counts = cn, samples = sam, study = list(id = "sim", feature_type = "otu", subject = "subj"))
  do.call(rbind, lapply(names(longc_engines), function(nm) {
    cfg <- longc_engines[[nm]]
    t0 <- proc.time()[["elapsed"]]
    e <- quiet_fit(st, ~ g + age, engine = cfg$engine, longitudinal = cfg$longitudinal, min_prevalence = 0, min_mean = 0,
                   exclude_unassigned = FALSE)
    secs <- proc.time()[["elapsed"]] - t0
    row <- if (is.null(e)) NULL else e[e$feature == "F1" & e$part == "abundance" & e$term == "gb", ]
    if (is.null(row) || !nrow(row)) return(data.frame(rep = r, engine = nm, estimate = NA, se = NA, p = NA, converged = NA, secs = secs))
    data.frame(rep = r, engine = nm, estimate = row$estimate, se = row$se, p = row$p_value, converged = row$converged, secs = secs)
  }))
}

grid <- switch(set, cross = cross_grid, long = long_grid, longc = longc_grid)
one <- switch(set, cross = cross_one, long = long_one, longc = longc_one)
only <- Sys.getenv("BENCH_ONLY")                                   # optional comma-separated scenario names
if (nzchar(only)) grid <- grid[grid$scenario %in% strsplit(only, ",")[[1]], , drop = FALSE]
for (i in seq_len(nrow(grid))) {
  sc <- grid[i, ]
  f <- file.path(out_dir, paste0(sc$scenario, ".csv"))
  if (file.exists(f)) next
  res <- do.call(rbind, mclapply(seq_len(nrep), one, sc = sc, mc.cores = cores))
  res <- cbind(scenario = sc$scenario, res, sc[setdiff(names(sc), "scenario")], row.names = NULL)
  utils::write.csv(res, f, row.names = FALSE)
  cat(sc$scenario, "done", format(Sys.time()), "\n")
}
writeLines(utils::capture.output(utils::sessionInfo()), file.path(out_dir, paste0("sessionInfo-", set, ".txt")))
cat("all done\n")
