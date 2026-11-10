# compare_mmuphin_melody.R
# Runs MMUPHin (lm_meta) and Melody (package miMeta) on the same 15 colorectal cancer
# studies that Claude analysed with metamicrobiomeR 2.0, so the three methods can be
# compared. Input: crc_data/crc_genus_tables.rds (genus counts with harmonised names,
# written by Claude). Output: crc_data/comparison_results.rds and comparison_log.txt.
#
# Run in a fresh R session:
#   setwd("C:/Users/ADMIN/Desktop/metamicrobiomeR/CRCstudies")
#   source("compare_mmuphin_melody.R")
# It takes about 10 to 20 minutes. Claude reads the results from the crc_data folder.

options(timeout = 1800)
dir <- "crc_data"
log_file <- file.path(dir, "comparison_log.txt")
zz <- file(log_file, open = "wt"); sink(zz, split = TRUE); sink(zz, type = "message")
on.exit({ sink(type = "message"); sink(); close(zz) }, add = TRUE)
cat("Started", format(Sys.time()), "\n")

## 1. Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("MMUPHin", quietly = TRUE)) BiocManager::install("MMUPHin", update = FALSE, ask = FALSE)
if (!requireNamespace("miMeta", quietly = TRUE)) remotes::install_github("ZjpWei/miMeta", upgrade = "never")
cat("MMUPHin", as.character(packageVersion("MMUPHin")), " miMeta",
    if (requireNamespace("miMeta", quietly = TRUE)) as.character(packageVersion("miMeta")) else "NOT INSTALLED", "\n")

## 2. One genus table for all studies
st <- readRDS(file.path(dir, "crc_genus_tables.rds"))
present <- table(unlist(lapply(st, function(o) rownames(o$counts))))
genera <- sort(names(present)[present >= length(st) / 2])     # in at least half of the studies, as in pool_effects()
cat("Studies:", length(st), " genera in at least half of them:", length(genera), "\n")
abd <- do.call(cbind, lapply(names(st), function(k) {
  m <- st[[k]]$counts; out <- matrix(0, length(genera), ncol(m), dimnames = list(genera, paste(k, colnames(m), sep = ":")))
  g <- intersect(genera, rownames(m)); out[g, ] <- m[g, ]; out }))
meta <- do.call(rbind, lapply(names(st), function(k) {
  s <- st[[k]]$samples; data.frame(sample_id = paste(k, rownames(s), sep = ":"), study = k, group = s$group,
                                     platform = st[[k]]$platform, stringsAsFactors = FALSE) }))
rownames(meta) <- meta$sample_id
meta$study <- factor(meta$study); meta$group <- factor(meta$group, levels = c("control", "CRC"))
rel <- sweep(abd, 2, pmax(colSums(abd), 1), "/")               # proportions over the kept genera
cat("Samples:", ncol(abd), "\n"); print(table(meta$study, meta$group))
res <- list(genera = genera, n_samples = ncol(abd))

## 3. MMUPHin, default settings (TSS, arcsine square root, linear models, REML pooling)
res$mmuphin <- tryCatch({
  f <- MMUPHin::lm_meta(feature_abd = rel, batch = "study", exposure = "group", data = meta,
                        control = list(verbose = FALSE, output = file.path(dir, "MMUPHin_lm_meta"), forest_plot = NULL))
  cat("MMUPHin done:", nrow(f$meta_fits), "genera\n"); f$meta_fits
}, error = function(e) { cat("MMUPHin ERROR:", conditionMessage(e), "\n"); conditionMessage(e) })

## 4. Melody (miMeta). The argument names are read from the installed package.
if (requireNamespace("miMeta", quietly = TRUE)) {
  fns <- ls(asNamespace("miMeta")); cat("miMeta functions:", paste(grep("melody", fns, value = TRUE), collapse = ", "), "\n")
  fm <- names(formals(miMeta::melody)); cat("melody() arguments:", paste(fm, collapse = ", "), "\n")
  res$melody_formals <- fm
  cand <- list(rel.abd = t(rel), sample.data = meta, sample.id = "sample_id", study = "study",
               disease = "group", covariate.interest = "group", verbose = TRUE)
  args <- cand[intersect(names(cand), fm)]
  try_melody <- function(a) tryCatch(do.call(miMeta::melody, a), error = function(e) e)
  out <- try_melody(args)
  if (inherits(out, "error")) {                                  # try the other orientation of the table
    cat("Melody with samples in rows failed:", conditionMessage(out), "\n")
    args$rel.abd <- rel; out <- try_melody(args)
  }
  if (inherits(out, "error")) { cat("Melody ERROR:", conditionMessage(out), "\n"); res$melody <- conditionMessage(out) }
  else { cat("Melody done\n"); res$melody <- out }
}

saveRDS(res, file.path(dir, "comparison_results.rds"))
cat("Finished", format(Sys.time()), "\nTell Claude it is done. Results: ", normalizePath(file.path(dir, "comparison_results.rds")), "\n")
