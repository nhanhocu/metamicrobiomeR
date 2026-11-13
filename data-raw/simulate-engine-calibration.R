# Usage: Rscript data-raw/simulate-engine-calibration.R <output.rds> <replicates>
# Results are in dev/engine-calibration-results.csv.
# Results are summarised in the description of the pull request that added bezi_glmmtmb and in the benchmark protocol.
# Null simulation through the package itself: fit_study(engine = "bezi_gamlss", longitudinal = TRUE) against
# fit_study(engine = "bezi_glmmtmb", longitudinal = TRUE). True group effect is zero throughout.
suppressPackageStartupMessages({library(metamicrobiomeR); library(parallel)})
options(metamicrobiomeR.quiet = TRUE)
a <- commandArgs(TRUE); out <- a[1]; nrep <- as.integer(a[2])
gr <- expand.grid(sdu = 0.5, prev = c(1, .6, .15), term = c("between", "within"), stringsAsFactors = FALSE)
sim <- function(de) {
  n <- 40; m <- 4; id <- rep(seq_len(n), each = m)
  g <- if (de$term == "between") rep(rep(c("a", "b"), length.out = n), each = m) else rep(c("a", "b"), length.out = n * m)
  age <- runif(n * m, 0, 6); u <- rnorm(n, sd = de$sdu)[id]
  mu <- plogis(-3 + u + .1 * age); phi <- 20
  y <- rbeta(n * m, mu * phi, (1 - mu) * phi); y[runif(n * m) > de$prev] <- 0
  sam <- data.frame(subj = paste0("s", id), g = g, age = age, row.names = paste0("S", seq_len(n * m)))
  rel <- rbind(F1 = y, other = 1 - y); colnames(rel) <- rownames(sam)
  as_study(relabund = rel, samples = sam, study = list(id = "sim", feature_type = "otu", subject = "subj"))
}
fit <- function(st, engine) {
  e <- suppressWarnings(try(fit_study(st, ~ g + age, engine = engine, longitudinal = TRUE, min_prevalence = 0, min_mean = 0,
                                      exclude_unassigned = FALSE), silent = TRUE))
  if (inherits(e, "try-error")) return(c(est = NA, se = NA, p = NA, conv = NA))
  r <- e[e$feature == "F1" & e$term == "gb" & e$part == "abundance", ]
  if (!nrow(r)) return(c(est = NA, se = NA, p = NA, conv = NA))
  c(est = r$estimate, se = r$se, p = r$p_value, conv = r$converged)
}
one <- function(r, de) {
  set.seed(5000 + r); st <- sim(de)
  rbind(data.frame(engine = "bezi_gamlss", t(fit(st, "bezi_gamlss"))), data.frame(engine = "bezi_glmmtmb", t(fit(st, "bezi_glmmtmb"))))
}
res <- list()
for (i in seq_len(nrow(gr))) {
  de <- gr[i, ]
  x <- do.call(rbind, mclapply(seq_len(nrep), one, de = de, mc.cores = 2))
  x$prev <- de$prev; x$term <- de$term; res[[i]] <- x
  saveRDS(do.call(rbind, res), out); cat(i, "of", nrow(gr), format(Sys.time()), "\n")
}
cat("done\n")
