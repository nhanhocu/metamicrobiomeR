S <- Sys.getenv("S"); base <- file.path(S, "crc", "hd_extract")
ranks <- c("kingdom","phylum","class","order","family","genus","species")
ids <- c(crc_baxter = "Baxter 2016", crc_zackular = "Zackular 2014", crc_zeller = "Zeller 2014", crc_zhao = "Wang 2012", crc_xiang = "Chen 2012")
out <- list()
for (id in names(ids)) {
  d <- file.path(base, paste0(id, "_results"), paste0(id, "_results"))
  otu <- read.delim(file.path(d, "RDP", paste0(id, ".otu_table.100.denovo.rdp_assigned")), row.names = 1, check.names = FALSE)
  l <- iconv(readLines(file.path(d, paste0(id, ".metadata.txt")), warn = FALSE), "latin1", "UTF-8")
  meta <- read.delim(text = l, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
  if (sum(rownames(otu) %in% rownames(meta)) > sum(colnames(otu) %in% rownames(meta))) otu <- t(otu)
  otu <- as.matrix(otu); sh <- intersect(colnames(otu), rownames(meta))
  cat(id, ": OTU table", nrow(otu), "x", ncol(otu), "| metadata", nrow(meta), "| matched", length(sh), "\n")
  print(table(meta[sh, "DiseaseState"], useNA = "ifany"))
  otu <- otu[, sh, drop = FALSE]; otu <- otu[rowSums(otu) > 0, , drop = FALSE]
  tax <- t(vapply(strsplit(rownames(otu), ";", fixed = TRUE), function(p) { v <- setNames(rep(NA_character_, 7), c("k","p","c","o","f","g","s"))
    for (q in p[grepl("^[kpcofgs]__", p)]) v[substr(q, 1, 1)] <- sub("^[a-z]__", "", q); v[v == ""] <- NA; v }, character(7)))
  tax <- as.data.frame(tax, stringsAsFactors = FALSE); names(tax) <- ranks
  out[[id]] <- list(counts = otu, taxonomy = tax, samples = meta[sh, , drop = FALSE], platform = "16S", label = ids[[id]])
  cat("  genera:", length(unique(na.omit(tax$genus))), " median reads:", median(colSums(otu)), " OTUs with genus:", round(mean(!is.na(tax$genus)) * 100), "%\n")
}
saveRDS(out, file.path(S, "crc", "amplicon_plain.rds"))
