# download_crc_studies.R
# Downloads public colorectal cancer (CRC) stool microbiome studies for the
# metamicrobiomeR 2.0 real-data example, and saves them in one file that Claude
# can read: crc_data/crc_studies.rds (plus a summary table).
#
# Sources
#   Shotgun (MetaPhlAn 3 species tables): curatedMetagenomicData (Bioconductor)
#   16S (OTU tables, RDP genus labels): MicrobiomeHD, Duvallet et al. 2017,
#     Zenodo record 840333 (CC BY-NC 4.0)
#
# How to run
#   1. Open a fresh R session (4.4 or newer), or in RStudio use Session > Restart R.
#   2. setwd() to a folder with about 1 GB free space.
#   3. source("download_crc_studies.R")
#   It takes 10 to 30 minutes, mostly for the first download from Bioconductor.
#   If it stops partway, run it again: finished steps are skipped.

options(timeout = 1800)
out_dir <- "crc_data"
dir.create(file.path(out_dir, "microbiomehd"), recursive = TRUE, showWarnings = FALSE)
log_line <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")

## 1. Packages -----------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("dplyr", "remotes")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

if (!requireNamespace("curatedMetagenomicData", quietly = TRUE)) {
  # curatedMetagenomicData 3.18.0 still calls rbiom::unifrac(), which rbiom 3.0 removed.
  # Install rbiom 2.2.1 (the last version with unifrac) first. It compiles from source,
  # so Windows needs Rtools, which you already have. Afterwards you can update rbiom
  # again with install.packages("rbiom") if another package needs the new version.
  rb <- if (requireNamespace("rbiom", quietly = TRUE)) utils::packageVersion("rbiom") else package_version("0")
  if (rb >= "3.0.0" || rb == "0") {
    log_line("Installing rbiom 2.2.1 (current rbiom is", as.character(rb), ")")
    if ("rbiom" %in% loadedNamespaces()) try(unloadNamespace("rbiom"), silent = TRUE)
    remotes::install_version("rbiom", version = "2.2.1", repos = "https://cloud.r-project.org",
                             upgrade = "never", dependencies = TRUE)
  }
  BiocManager::install(c("SummarizedExperiment", "TreeSummarizedExperiment", "curatedMetagenomicData"),
                       update = FALSE, ask = FALSE)
}
if (!requireNamespace("curatedMetagenomicData", quietly = TRUE))
  stop("curatedMetagenomicData could not be installed. Please send Claude the messages above.")
suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
})

## 2. Shotgun studies from curatedMetagenomicData ---------------------------------
shotgun_file <- file.path(out_dir, "shotgun_raw.rds")
if (!file.exists(shotgun_file)) {
  log_line("Finding CRC studies in curatedMetagenomicData")
  md <- curatedMetagenomicData::sampleMetadata
  md <- md[md$body_site == "stool" & md$study_condition %in% c("CRC", "control"), ]
  # keep studies with at least 10 cases and 10 controls
  tab <- table(md$study_name, md$study_condition)
  keep <- rownames(tab)[tab[, "CRC"] >= 10 & tab[, "control"] >= 10]
  md <- md[md$study_name %in% keep, ]
  md <- md[, colSums(!is.na(md)) > 0]
  log_line("Studies:", paste(keep, collapse = ", "))
  log_line("Downloading species tables for", nrow(md), "samples")
  tse <- curatedMetagenomicData::returnSamples(md, "relative_abundance", counts = TRUE, rownames = "long")
  saveRDS(tse, shotgun_file)
}
tse <- readRDS(shotgun_file)

ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
split_long <- function(x) {
  # "k__Bacteria|p__Firmicutes|...|s__Faecalibacterium_prausnitzii"
  parts <- strsplit(x, "|", fixed = TRUE)
  t(vapply(parts, function(p) {
    p <- sub("^[a-z]__", "", p)
    c(p, rep(NA_character_, 7))[1:7]
  }, character(7)))
}

shotgun <- list()
cd <- as.data.frame(colData(tse))
for (s in sort(unique(cd$study_name))) {
  idx <- which(cd$study_name == s)
  cnt <- as.matrix(assay(tse)[, idx, drop = FALSE])
  cnt <- cnt[rowSums(cnt) > 0, , drop = FALSE]
  tax <- as.data.frame(split_long(rownames(cnt)), stringsAsFactors = FALSE)
  names(tax) <- ranks
  samp <- cd[idx, , drop = FALSE]
  samp <- samp[, colSums(!is.na(samp)) > 0, drop = FALSE]
  samp$sample_id <- colnames(cnt)
  shotgun[[s]] <- list(
    counts = cnt, taxonomy = tax, samples = samp,
    platform = "shotgun", profiler = "MetaPhlAn 3 (curatedMetagenomicData)",
    country = paste(unique(samp$country), collapse = ";"),
    source = "curatedMetagenomicData"
  )
}
log_line("Shotgun studies prepared:", length(shotgun))

## 3. 16S studies from MicrobiomeHD -----------------------------------------------
hd_files <- c(crc_baxter = "crc_baxter_results.tar.gz", crc_zackular = "crc_zackular_results.tar.gz",
              crc_zeller = "crc_zeller_results.tar.gz", crc_wang = "crc_zhao_results.tar.gz",
              crc_chen = "crc_xiang_results.tar.gz", crc_zhu = "crc_zhu_results.tar.gz")
amplicon <- list()
for (id in names(hd_files)) {
  f <- file.path(out_dir, "microbiomehd", hd_files[[id]])
  if (!file.exists(f) || file.size(f) < 100) {
    log_line("Downloading", hd_files[[id]])
    url <- paste0("https://zenodo.org/records/840333/files/", hd_files[[id]], "?download=1")
    ok <- tryCatch({ utils::download.file(url, f, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) { log_line("  failed:", conditionMessage(e)); FALSE })
    if (!ok) next
  }
  ex <- file.path(out_dir, "microbiomehd", id)
  if (!dir.exists(ex)) utils::untar(f, exdir = ex)
  all <- list.files(ex, recursive = TRUE, full.names = TRUE)
  otu_f <- grep("otu_table.*rdp_assigned$", all, value = TRUE)
  meta_f <- grep("metadata\\.txt$", all, value = TRUE)
  if (!length(otu_f) || !length(meta_f)) {
    log_line("  could not find the OTU table or metadata for", id, "- files are kept for Claude to inspect")
    amplicon[[id]] <- list(files = sub(paste0("^", ex, "/"), "", all))
    next
  }
  otu <- utils::read.delim(otu_f[1], row.names = 1, check.names = FALSE)
  meta <- utils::read.delim(meta_f[1], row.names = 1, check.names = FALSE,
                            fileEncoding = "latin1", stringsAsFactors = FALSE)
  # orient as taxa x samples, matching the metadata sample IDs
  if (sum(rownames(otu) %in% rownames(meta)) > sum(colnames(otu) %in% rownames(meta))) otu <- t(otu)
  otu <- as.matrix(otu)
  shared <- intersect(colnames(otu), rownames(meta))
  otu <- otu[, shared, drop = FALSE]
  otu <- otu[rowSums(otu) > 0, , drop = FALSE]
  # RDP labels look like k__Bacteria;p__Firmicutes;...;g__Blautia;s__;d__denovo123
  lab <- strsplit(rownames(otu), ";", fixed = TRUE)
  tax <- t(vapply(lab, function(p) {
    p <- p[grepl("^[kpcofgs]__", p)]
    v <- setNames(rep(NA_character_, 7), c("k", "p", "c", "o", "f", "g", "s"))
    for (q in p) v[substr(q, 1, 1)] <- sub("^[a-z]__", "", q)
    v[v == ""] <- NA
    v
  }, character(7)))
  tax <- as.data.frame(tax, stringsAsFactors = FALSE)
  names(tax) <- ranks
  samp <- meta[shared, , drop = FALSE]
  samp$sample_id <- shared
  amplicon[[id]] <- list(
    counts = otu, taxonomy = tax, samples = samp,
    platform = "16S", profiler = "MicrobiomeHD pipeline, RDP genus labels",
    source = "MicrobiomeHD (Zenodo 840333)", raw_files = basename(c(otu_f[1], meta_f[1]))
  )
  log_line("  ", id, ":", ncol(otu), "samples,", nrow(otu), "OTUs")
}

## 4. Save one file and a summary ----------------------------------------------------
crc <- list(shotgun = shotgun, amplicon = amplicon,
            created = Sys.time(), session = utils::sessionInfo())
saveRDS(crc, file.path(out_dir, "crc_studies.rds"), compress = "xz")

`%||%` <- function(a, b) if (is.null(a)) b else a
describe <- function(x, name, type) {
  if (is.null(x$counts)) return(data.frame(study = name, type = type, samples = NA, taxa = NA,
                                           groups = "not parsed", median_depth = NA))
  g <- x$samples$study_condition %||% x$samples$DiseaseState
  data.frame(study = name, type = type, samples = ncol(x$counts), taxa = nrow(x$counts),
             groups = if (is.null(g)) "unknown" else paste(names(table(g)), table(g), sep = " ", collapse = ", "),
             median_depth = stats::median(colSums(x$counts)))
}
summ <- rbind(do.call(rbind, Map(describe, shotgun, names(shotgun), "shotgun")),
              do.call(rbind, Map(describe, amplicon, names(amplicon), "16S")))
utils::write.csv(summ, file.path(out_dir, "crc_studies_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)
log_line("Done. Please send Claude these two files:")
log_line("  ", normalizePath(file.path(out_dir, "crc_studies.rds")),
         sprintf("(%.1f MB)", file.size(file.path(out_dir, "crc_studies.rds")) / 1e6))
log_line("  ", normalizePath(file.path(out_dir, "crc_studies_summary.csv")))
