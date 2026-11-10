# download_16s_microbiomehd.R
# Downloads the six colorectal cancer 16S studies of MicrobiomeHD (Duvallet et al. 2017,
# Zenodo record 840333, CC BY-NC 4.0) into crc_data/microbiomehd. Claude reads and
# processes the archives, so nothing else is needed here.
#
# Run it in the same folder as before (the one that holds crc_data):
#   setwd("C:/Users/ADMIN/Desktop/metamicrobiomeR/CRCstudies")
#   source("download_16s_microbiomehd.R")

options(timeout = 1800)
dest <- file.path("crc_data", "microbiomehd")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

files <- c("crc_baxter_results.tar.gz",   # Baxter 2016, USA, V4, about 18 MB
           "crc_zeller_results.tar.gz",   # Zeller 2014, France, V4, about 12 MB
           "crc_zackular_results.tar.gz", # Zackular 2014, USA, V4, about 8 MB
           "crc_zhu_results.tar.gz",      # Zhu, China, under 1 MB
           "crc_zhao_results.tar.gz",     # Wang 2012, China, under 1 MB
           "crc_xiang_results.tar.gz")    # Chen 2012, China, under 1 MB

for (f in files) {
  out <- file.path(dest, f)
  if (file.exists(out) && file.size(out) > 1000) {
    cat("already here:", f, "\n"); next
  }
  url <- paste0("https://zenodo.org/records/840333/files/", f, "?download=1")
  cat("downloading", f, "... ")
  res <- tryCatch(utils::download.file(url, out, mode = "wb", quiet = TRUE, method = "libcurl"),
                  error = function(e) conditionMessage(e), warning = function(w) conditionMessage(w))
  if (identical(res, 0L) && file.exists(out)) {
    cat(sprintf("done (%.1f MB)\n", file.size(out) / 1e6))
  } else {
    cat("FAILED:", as.character(res), "\n")
    if (file.exists(out) && file.size(out) < 1000) file.remove(out)
  }
}

cat("\nFiles in", normalizePath(dest), ":\n")
print(data.frame(file = list.files(dest), MB = round(file.size(list.files(dest, full.names = TRUE)) / 1e6, 2)))
cat("\nIf any file FAILED, open this page in your browser and download it into the folder above:\n",
    "https://zenodo.org/records/840333\n")
