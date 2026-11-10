suppressMessages(library(metamicrobiomeR))
S <- Sys.getenv("S"); st <- readRDS(file.path(S, "crc", "crc_studies_genus.rds"))
main <- st[names(st) != "ZellerG_2014_16S"]
runs <- list(
  list(name = "clr",    engine = "lm_clr",       f = ~ group,             set = names(main)),
  list(name = "asin",   engine = "lm_asin",      f = ~ group,             set = names(main)),
  list(name = "zibb",   engine = "zibb_glmmtmb", f = ~ group,             set = names(main)),
  list(name = "clr_adj", engine = "lm_clr",      f = ~ group + age + sex, set = setdiff(names(main), c("WangT_2012_16S", "ChenW_2012_16S"))),
  list(name = "bezi",   engine = "bezi_gamlss",  f = ~ group,             set = names(main)),
  list(name = "zibb_adj", engine = "zibb_glmmtmb", f = ~ group + age + sex, set = setdiff(names(main), c("WangT_2012_16S", "ChenW_2012_16S"))))
for (r in runs) for (k in r$set) {
  out <- file.path(S, "crc", "fits", paste0(r$name, "__", k, ".rds"))
  if (file.exists(out)) next
  t0 <- Sys.time(); w <- character()
  x <- st[[k]]
  if (grepl("age", deparse(r$f))) x <- subset_samples(x, stats::complete.cases(mm_samples(x)[, c("age", "sex")]))
  e <- withCallingHandlers(tryCatch(fit_study(x, r$f, engine = r$engine, cores = 2), error = function(err) { message("ERROR ", k, " ", r$name, ": ", conditionMessage(err)); NULL }),
                           warning = function(wn) { w <<- c(w, conditionMessage(wn)); invokeRestart("muffleWarning") }, message = function(m) invokeRestart("muffleMessage"))
  attr(e, "warnings_captured") <- unique(w)
  saveRDS(e, out)
  cat(format(Sys.time(), "%H:%M:%S"), r$name, k, if (is.null(e)) "FAILED" else nrow(e), "rows,", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}
cat("all fits done\n")
