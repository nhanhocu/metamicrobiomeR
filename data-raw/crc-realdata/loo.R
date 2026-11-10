suppressMessages(library(metamicrobiomeR))
S <- Sys.getenv("S"); fd <- file.path(S, "crc", "fits"); st <- readRDS(file.path(S, "crc", "crc_studies_genus.rds"))
ids <- setdiff(names(st), "ZellerG_2014_16S")
q <- function(expr) withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"), message = function(m) invokeRestart("muffleMessage"))
for (e in c("clr", "bezitmb")) {
  eff <- do.call(bind_effects, lapply(ids, function(k) readRDS(file.path(fd, paste0(e, "__", k, ".rds")))))
  eff <- structure(q(filter_effects(eff)), class = class(eff))
  full <- as.data.frame(q(pool_effects(eff, term = "groupCRC")))
  sig <- full[full$p_adjust < 0.05, c("feature", "part")]
  l <- as.data.frame(q(leave_one_out(eff, term = "groupCRC")))
  key <- paste(sig$feature, sig$part)
  l <- l[paste(l$feature, l$part) %in% key, ]
  res <- do.call(rbind, lapply(split(l, paste(l$feature, l$part)), function(d) data.frame(key = paste(d$feature[1], d$part[1]), omissions = nrow(d), still_sig = sum(d$p_adjust < 0.05), worst_p = max(d$p_adjust), worst_without = d$omitted[which.max(d$p_adjust)])))
  cat("==", e, ":", nrow(sig), "significant;", sum(res$still_sig == res$omissions), "stay significant whichever study is left out\n")
  print(res[res$still_sig < res$omissions, ], row.names = FALSE)
  saveRDS(res, file.path(S, "crc", "pooled", paste0("loo_", e, ".rds")))
}
