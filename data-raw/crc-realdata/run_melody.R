# run_melody.R
# Runs Melody (package miMeta) on the same 15 colorectal cancer studies.
# The installed melody() takes feature.counts (one table per study) and
# covariate.interest (one data frame per study). This script also saves the help
# page and the structure of the result, so Claude can check how it was used.
#
#   setwd("C:/Users/ADMIN/Desktop/metamicrobiomeR/CRCstudies")
#   source("crc_data/run_melody.R")

dir <- "crc_data"
log_file <- file.path(dir, "melody_log.txt")
zz <- file(log_file, open = "wt"); sink(zz, split = TRUE)
on.exit(sink(), add = TRUE); on.exit(close(zz), add = TRUE)
library(miMeta)
cat("miMeta", as.character(packageVersion("miMeta")), "\n")

# 1. Save the help page of melody() as text, for Claude
hf <- utils:::.getHelpFile(utils::help("melody", package = "miMeta"))
tools::Rd2txt(hf, out = file.path(dir, "melody_help.txt"), options = list(underline_titles = FALSE))
cat("Help page saved\n")

# 2. Inputs: one sample-by-genus count table and one data frame per study
st <- readRDS(file.path(dir, "crc_genus_tables.rds"))
present <- table(unlist(lapply(st, function(o) rownames(o$counts))))
genera <- sort(names(present)[present >= length(st) / 2])   # the same 158 genera as for MMUPHin
feature.counts <- list(); covariate.interest <- list()
for (k in names(st)) {
  m <- st[[k]]$counts
  out <- matrix(0, ncol(m), length(genera), dimnames = list(paste(k, colnames(m), sep = ":"), genera))
  g <- intersect(genera, rownames(m)); out[, g] <- t(m[g, , drop = FALSE])
  feature.counts[[k]] <- out
  covariate.interest[[k]] <- data.frame(disease = as.numeric(st[[k]]$samples$group == "CRC"), row.names = rownames(out))
}
cat("Studies:", length(feature.counts), " genera:", length(genera), " samples:", sum(sapply(feature.counts, nrow)), "\n")

# 3. Melody with its default settings
t0 <- Sys.time()
res <- tryCatch(melody(feature.counts = feature.counts, covariate.interest = covariate.interest, verbose = TRUE),
                error = function(e) e)
cat("Time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutes\n")
if (inherits(res, "error")) {
  cat("Melody ERROR:", conditionMessage(res), "\n")
} else {
  cat("Melody done. Structure of the result:\n"); str(res, max.level = 3)
}
saveRDS(list(result = res, genera = genera), file.path(dir, "melody_results.rds"))
cat("Saved", normalizePath(file.path(dir, "melody_results.rds")), "\nTell Claude it is done.\n")
