suppressMessages(library(metamicrobiomeR))
S <- Sys.getenv("S"); st <- readRDS(file.path(S, "crc", "crc_studies_genus.rds")); main <- st[names(st) != "ZellerG_2014_16S"]
for (k in names(main)) {
  out <- file.path(S, "crc", "fits", paste0("bezitmb__", k, ".rds")); if (file.exists(out)) next
  t0 <- Sys.time(); w <- character()
  e <- withCallingHandlers(tryCatch(fit_study(main[[k]], ~ group, engine = "bezi_glmmtmb", prevalence = ~ group, cores = 2), error = function(err) { message("ERROR ", k, ": ", conditionMessage(err)); NULL }),
         warning = function(wn) { w <<- c(w, conditionMessage(wn)); invokeRestart("muffleWarning") }, message = function(m) invokeRestart("muffleMessage"))
  attr(e, "warnings_captured") <- unique(w); saveRDS(e, out)
  x <- as.data.frame(e); cat(format(Sys.time(), "%H:%M:%S"), k, nrow(x), "rows,", sum(x$converged), "converged,", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")
}
cat("bezitmb done\n")
