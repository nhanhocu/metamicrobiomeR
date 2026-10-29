# Diagnostics for the unexplained standard-error differences in the stored Bangladesh sex results
# (see dev/VALIDATION.md, F11). oldfit.R needs gamlss 5.0-5 and gamlss.dist 5.0-3 installed in a library
# that comes first in R_LIBS_USER; qrtest.R compares summary(type = 'vcov') with type = 'qr'.
# Both expect dev/validation/bangladesh-sex-fit.rds from data-raw/refit-sex-bangladesh.R.
suppressMessages(library(gamlss)); source("data-raw/legacy-helpers.R")
env <- new.env(); for (f in list.files("data-raw/legacy-1.1/results", full.names = TRUE)) load(f, env)
s <- legacy_to_study(as.data.frame(get("taxtab6.rm", env)[[5]]), list(id = "B", feature_type = "otu", subject = "personid"))$study
st <- get("taxacom6.zi.rm.sex.adjustbfage", env)
sd <- mm_samples(s); ra <- mm_relabund(s)
new <- readRDS("dev/validation/bangladesh-sex-fit.rds"); n <- new[new$term == "genderMale", ]; n <- n[match(st$id, n$feature), ]
rel <- abs(n$se - st[["Std. Error.genderMale"]]) / st[["Std. Error.genderMale"]]
set.seed(2); ids <- c(st$id[order(-rel)][1:6], sample(st$id[rel <= 1e-6], 4))
out <- lapply(ids, function(id) {
  d <- data.frame(.y = ra[id, ], gender = factor(sd$gender), bf = factor(sd$bf), age.sample = sd$age.sample, .subject = factor(sd$personid))
  d$.y[which(d$.y >= 1)] <- 0.9999; d <- na.omit(d)
  fit <- suppressWarnings(gamlss(.y ~ gender + bf + age.sample + random(.subject), family = "BEZI", data = d, trace = FALSE, control = gamlss.control(n.cyc = 20, trace = FALSE)))
  get_se <- function(type) { msgs <- character(); x <- NULL
    capture.output(x <- withCallingHandlers(try(summary(fit, type = type, save = TRUE), silent = TRUE), warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") }))
    list(se = if (inherits(x, "try-error")) NA else x$coef.table["genderMale", "Std. Error"], vcov_failed = any(grepl("vcov has failed", msgs))) }
  v <- get_se("vcov"); q <- get_se("qr")
  data.frame(se_stored = st[["Std. Error.genderMale"]][st$id == id], se_vcov = v$se, se_qr = q$se, vcov_fell_back_to_qr = v$vcov_failed)
})
r <- cbind(group = rep(c("worst", "matching"), c(6, 4)), do.call(rbind, out)); print(r, digits = 4, row.names = FALSE)
cat(sprintf("max |qr - stored| = %.1e ; max |vcov - stored| = %.1e\n", max(abs(r$se_qr - r$se_stored)), max(abs(r$se_vcov - r$se_stored))))
