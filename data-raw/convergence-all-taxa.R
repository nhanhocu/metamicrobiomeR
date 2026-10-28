# Example 1, all ranks: how much do estimates change when non-converged GAMLSS fits are allowed to finish?
source("data-raw/legacy-helpers.R")
env <- new.env(); for (f in list.files("data-raw/legacy-1.1/results", full.names = TRUE)) load(f, env)
s <- legacy_to_study(as.data.frame(get("taxtab6.rm", env)[[5]]), list(id = "ex1", feature_type = "otu", subject = "personid"))$study
f <- ~ bf + age.sample
base <- fit_study(s, f, longitudinal = TRUE, n_cycles = 20, exclude_unassigned = FALSE)
info <- attr(base, "fit_info"); bad <- info$feature[!info$converged]
set.seed(1); ids <- sample(bad, 25)
sub <- s; sub$assays <- lapply(s$assays, function(m) m[ids, , drop = FALSE]); sub$features <- s$features[ids, , drop = FALSE]
long <- fit_study(sub, f, longitudinal = TRUE, n_cycles = 300, exclude_unassigned = FALSE)
a <- base[base$feature %in% ids, ]; m <- merge(a, long, by = c("feature", "term"), suffixes = c("_20", "_300"))
m$d_est <- m$estimate_300 - m$estimate_20; m$d_p <- m$p_value_300 - m$p_value_20
m$sig_changes <- (m$p_value_20 < 0.05) != (m$p_value_300 < 0.05)
li <- attr(long, "fit_info")
cat(sprintf("Non-converged at 20 cycles: %d of %d features; sampled %d\n", length(bad), nrow(info), length(ids)))
cat(sprintf("  of the sample, converged with 300 cycles: %d of %d\n", sum(li$converged), nrow(li)))
for (tm in unique(m$term)) { x <- m[m$term == tm, ]
  cat(sprintf("  %-18s median |change| estimate %.4f (max %.3f); median |change| p %.4f (max %.3f); significance flips: %d of %d\n",
      tm, median(abs(x$d_est), na.rm = TRUE), max(abs(x$d_est), na.rm = TRUE), median(abs(x$d_p), na.rm = TRUE), max(abs(x$d_p), na.rm = TRUE), sum(x$sig_changes, na.rm = TRUE), nrow(x))) }
