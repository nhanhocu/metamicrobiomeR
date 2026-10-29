# Summarise the pilot benchmark (data-raw/benchmark-pilot.R). Run from the repository root.
fs <- list.files("dev/validation/bench", pattern = "\\.rds$", full.names = TRUE)
d <- do.call(rbind, lapply(fs, readRDS))
summ <- function(x) {
  x$sig <- !is.na(x$p) & x$p < 0.05
  # BH within each replicate (one study, all features), target 5%
  x$bh <- ave(x$p, x$rep, FUN = function(p) p.adjust(p, "fdr")) < 0.05
  nd <- x[!x$diff, ]; df <- x[x$diff, ]
  rep_fdp <- tapply(seq_len(nrow(x)), x$rep, function(i) { r <- i[x$bh[i]]; if (length(r)) mean(!x$diff[r]) else NA })
  data.frame(
    reps = length(unique(x$rep)),
    false_positive_rate = if (nrow(nd)) mean(nd$sig, na.rm = TRUE) else NA,
    power = if (nrow(df)) mean(df$sig, na.rm = TRUE) else NA,
    mean_FDP_at_BH5 = if (any(x$diff)) mean(rep_fdp, na.rm = TRUE) else NA,
    converged = mean(x$converged),
    sec_per_rep = mean(tapply(x$seconds, x$rep, mean)))
}
key <- paste(d$scenario, d$engine, sep = "|")
out <- do.call(rbind, lapply(split(d, key), summ))
out <- cbind(scenario = sub("\\|.*", "", rownames(out)), engine = sub(".*\\|", "", rownames(out)), out)
out <- out[order(match(out$scenario, c("null", "alt_balanced", "alt_unbalanced")), out$engine), ]; rownames(out) <- NULL
num <- vapply(out, is.numeric, logical(1)); out[num] <- lapply(out[num], function(v) round(v, 3))
print(out, row.names = FALSE)
# Monte Carlo uncertainty of a rate: binomial standard error
n_null <- tapply(!d$diff[d$scenario == "null"], d$engine[d$scenario == "null"], sum)
cat("\nNull scenario: tests per engine:", paste(names(n_null), n_null, collapse = "; "), "\n")
cat("Approx. 95% Monte Carlo interval for a 5% rate with that many tests: 0.05 +/-", round(1.96 * sqrt(0.05 * 0.95 / min(n_null)), 3), "\n")
