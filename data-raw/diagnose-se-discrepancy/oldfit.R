# Diagnostics for the unexplained standard-error differences in the stored Bangladesh sex results
# (see dev/VALIDATION.md, F11). oldfit.R needs gamlss 5.0-5 and gamlss.dist 5.0-3 installed in a library
# that comes first in R_LIBS_USER; qrtest.R compares summary(type = 'vcov') with type = 'qr'.
# Both expect dev/validation/bangladesh-sex-fit.rds from data-raw/refit-sex-bangladesh.R.
source("data-raw/legacy-helpers.R")
cat("gamlss in use:", as.character(packageVersion("gamlss")), "\n")
env <- new.env(); for (f in list.files("data-raw/legacy-1.1/results", full.names = TRUE)) load(f, env)
s <- legacy_to_study(as.data.frame(get("taxtab6.rm", env)[[5]]), list(id = "Bangladesh", feature_type = "otu", subject = "personid"))$study
new <- readRDS("dev/validation/bangladesh-sex-fit.rds"); st <- get("taxacom6.zi.rm.sex.adjustbfage", env)
n <- new[new$term == "genderMale", ]; n <- n[match(st$id, n$feature), ]
rel <- abs(n$se - st[["Std. Error.genderMale"]]) / st[["Std. Error.genderMale"]]
set.seed(2); ids <- c(st$id[order(-rel)][1:6], sample(st$id[rel <= 1e-6], 4))
sub <- s; sub$assays <- lapply(s$assays, function(m) m[ids, , drop = FALSE]); sub$features <- s$features[ids, , drop = FALSE]
t0 <- Sys.time()
old <- suppressWarnings(fit_study(sub, ~ gender + bf + age.sample, longitudinal = TRUE, n_cycles = 20, exclude_unassigned = FALSE))
cat(sprintf("refit %d taxa in %.0f s\n", length(ids), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
o <- old[old$term == "genderMale", ]; o <- o[match(ids, o$feature), ]
cmp <- data.frame(group = rep(c("worst", "matching"), c(6, 4)), se_stored = st[["Std. Error.genderMale"]][match(ids, st$id)],
                  se_gamlss_5.5 = n$se[match(ids, n$feature)], se_gamlss_5.0.5 = o$se,
                  est_diff_old_vs_stored = abs(o$estimate - st[["Estimate.genderMale"]][match(ids, st$id)]))
print(cmp, digits = 4, row.names = FALSE)
cat(sprintf("max |se(5.0-5) - se(stored)| = %.1e ; max |se(5.5-5) - se(stored)| = %.1e\n", max(abs(cmp$se_gamlss_5.0.5 - cmp$se_stored)), max(abs(cmp$se_gamlss_5.5 - cmp$se_stored))))
