# (n_cycles = 20, the 1.x default, is used throughout to reproduce 1.x.)
# Reproduce the GAMLSS-BEZI results of Ho et al. (2019), Tables 3 to 5, with the 2.0 engine.
# Compares (a) the new engine with the stored 1.x results objects and (b) with the
# numbers printed in the paper. Run from the repository root:
#   Rscript data-raw/validate-2019-tables.R          # phylum level (fast)
#   FULL=1 Rscript data-raw/validate-2019-tables.R   # also all taxa of Example 1
source("data-raw/legacy-helpers.R")
env <- new.env()
for (f in list.files("data-raw/legacy-1.1/results", full.names = TRUE)) load(f, env)

phyla <- c("actinobacteria", "bacteroidetes", "firmicutes", "proteobacteria")
P <- function(ex, grp, est, lo, hi, p, fdr) data.frame(example = ex, phylum = phyla, est, lo, hi, p, fdr)
paper <- rbind(
  P("1",  "non-EBF vs EBF",       c(-0.37, 0.26, 0.24, 0.37), c(-0.65, 0.00, 0.00, 0.11), c(-0.10, 0.53, 0.47, 0.64),
    c(.0083, .0499, .0468, .0053), c(.0166, .0499, .0499, .0166)),
  P("2",  "solid food >5 months", c(0.19, -0.26, -0.16, 0.14), c(0.04, -0.42, -0.30, -0.02), c(0.34, -0.10, -0.03, 0.30),
    c(.0119, .0018, .0156, .0861), c(.0208, .0070, .0208, .0861)),
  P("3a", "diarrhea, EBF <=2 mo", c(-0.73, -0.29, 0.49, -0.17), c(-1.12, -0.68, 0.15, -0.54), c(-0.34, 0.10, 0.84, 0.20),
    c(.0003, .1524, .0055, .3729), c(.0011, .2032, .0109, .3729)),
  P("3b", "diarrhea, EBF >2 mo",  c(0.02, 0.07, -0.02, 0.12), c(-0.42, -0.41, -0.40, -0.33), c(0.46, 0.56, 0.36, 0.56),
    c(.9243, .7680, .9142, .6043), c(.9243, .9243, .9243, .9243)))

examples <- list(
  "1"  = list(tab = "taxtab6.rm",            formula = ~ bf + age.sample,          term = "bfNon_exclusiveBF",    long = TRUE,
              stored = "taxacom6.zi.rmg"),
  "2"  = list(tab = "taxtab6plus.rm",        formula = ~ month.food5 + age.sample, term = "month.food5>5 months", long = TRUE,
              stored = "taxacom.6plus.sl5.zi.rmg"),
  # The stored results of Example 3 were computed with longitudinal = "no" (see the
  # supplementary file), although the paper describes longitudinal data.
  "3a" = list(tab = "taxtab6plus.exbf2.rm",     formula = ~ diarrhea + age.sample, term = "diarrheaYes", long = FALSE,
              stored = "taxacom.6plus.dia.exbf2.zi.rmg"),
  "3b" = list(tab = "taxtab6plus.exbf2plus.rm", formula = ~ diarrhea + age.sample, term = "diarrheaYes", long = FALSE,
              stored = "taxacom.6plus.dia.exbf2plus.zi.rmg"))

run <- function(ex, rank = "phylum") {
  e <- examples[[ex]]
  tab <- as.data.frame(get(e$tab, env)[[5]])
  lt <- legacy_to_study(tab, list(id = paste0("bangladesh_ex", ex), feature_type = "otu", subject = "personid"))
  fit <- fit_study(lt$study, e$formula, rank = rank, longitudinal = e$long, n_cycles = 20, exclude_unassigned = FALSE)
  fit <- fit[fit$term == e$term, ]
  fit$lo <- fit$estimate - 1.96 * fit$se; fit$hi <- fit$estimate + 1.96 * fit$se
  list(fit = fit, stored = get(e$stored, env), info = attr(fit, "fit_info"))
}

short <- function(x) sub("^k__bacteria\\.p__", "", x)
rows <- list(); diffs <- list()
for (ex in names(examples)) {
  t0 <- Sys.time(); r <- run(ex)
  f <- r$fit; f$phylum <- short(f$feature)
  f <- f[f$phylum %in% phyla, ]
  f$fdr <- stats::p.adjust(f$p_value, "fdr")      # BH within the four main phyla, as displayed in the paper
  # (a) against the stored 1.x results
  st <- r$stored; st$phylum <- short(st$id)
  tm <- examples[[ex]]$term
  m <- merge(f, st, by = "phylum")
  diffs[[ex]] <- data.frame(example = ex, n_phyla = nrow(m),
    max_abs_diff_estimate = max(abs(m$estimate - m[[paste0("Estimate.", tm)]])),
    max_abs_diff_se = max(abs(m$se - m[[paste0("Std. Error.", tm)]])),
    max_abs_diff_p = max(abs(m$p_value - m[[paste0("Pr(>|t|).", tm)]])))
  # (b) against the paper (values printed to 2 decimals, p-values to 4)
  pp <- paper[paper$example == ex, ]
  names(pp)[3:7] <- paste0("paper_", names(pp)[3:7])
  nw <- f[, c("phylum", "estimate", "lo", "hi", "p_value", "fdr")]
  names(nw)[2:6] <- c("new_est", "new_lo", "new_hi", "new_p", "new_fdr")
  g <- merge(pp, nw, by = "phylum")
  rows[[ex]] <- data.frame(example = ex, phylum = g$phylum,
    est = abs(round(g$new_est, 2) - g$paper_est) < 0.0051, lo = abs(round(g$new_lo, 2) - g$paper_lo) < 0.0051,
    hi = abs(round(g$new_hi, 2) - g$paper_hi) < 0.0051, p = abs(round(g$new_p, 4) - g$paper_p) < 0.00005,
    fdr = abs(round(g$new_fdr, 4) - g$paper_fdr) < 0.00005,
    new_est = round(g$new_est, 2), paper_est = g$paper_est, new_p = round(g$new_p, 4), paper_p = g$paper_p,
    new_fdr = round(g$new_fdr, 4), paper_fdr = g$paper_fdr)
  cat(sprintf("Example %s done in %.0f s; converged %d/%d features\n", ex, as.numeric(difftime(Sys.time(), t0, units = "secs")),
              sum(r$info$converged), nrow(r$info)))
}
cat("\n== (a) new engine vs stored 1.x results (phylum level) ==\n")
print(do.call(rbind, diffs), row.names = FALSE, digits = 3)
res <- do.call(rbind, rows)
cat("\n== (b) new engine vs numbers printed in the paper (Tables 3-5, GAMLSS-BEZI) ==\n")
print(res, row.names = FALSE)
cat(sprintf("\nCells matching the paper (estimate, CI limits, p, FDR): %d of %d\n",
            sum(as.matrix(res[, c("est", "lo", "hi", "p", "fdr")])), 5 * nrow(res)))
invisible(list(diffs = diffs, res = res))

# Optional: every taxon of Example 1 (all ranks), against the stored 1.x results
if (Sys.getenv("FULL") == "1") {
  t0 <- Sys.time()
  e <- examples[["1"]]
  tab <- as.data.frame(get(e$tab, env)[[5]])
  s <- legacy_to_study(tab, list(id = "bangladesh_ex1", feature_type = "otu", subject = "personid"))$study
  fit <- fit_study(s, e$formula, longitudinal = e$long, n_cycles = 20, exclude_unassigned = FALSE)
  st <- get(e$stored, env)
  cat(sprintf("\n== Example 1, all taxa: %d fitted by 2.0, %d in the stored 1.x result (%.0f s) ==\n",
              length(unique(fit$feature)), nrow(st), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  cat("Same set of taxa:", setequal(unique(fit$feature), st$id), "\n")
  for (tm in c("bfNon_exclusiveBF", "bfNo_BF", "age.sample")) {
    x <- fit[fit$term == tm, ]; x <- x[match(st$id, x$feature), ]
    cat(sprintf("  %-20s max |diff| estimate %.1e, se %.1e, p %.1e\n", tm,
                max(abs(x$estimate - st[[paste0("Estimate.", tm)]]), na.rm = TRUE),
                max(abs(x$se - st[[paste0("Std. Error.", tm)]]), na.rm = TRUE),
                max(abs(x$p_value - st[[paste0("Pr(>|t|).", tm)]]), na.rm = TRUE)))
  }
  info <- attr(fit, "fit_info")
  cat(sprintf("  converged at 20 cycles: %d of %d features\n", sum(info$converged), nrow(info)))
}
