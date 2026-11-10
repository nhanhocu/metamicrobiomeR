x <- readRDS("/mnt/user-data/uploads/crc_data/shotgun_raw.rds")
a <- attributes(x); cd <- attributes(a$colData)
samp <- as.data.frame(cd$listData, stringsAsFactors = FALSE); rownames(samp) <- cd$rownames
m <- attributes(attributes(a$assays)$data)$listData[[1]]
stopifnot(identical(colnames(m), rownames(samp)))
ranks <- c("kingdom","phylum","class","order","family","genus","species")
out <- list()
for (s in sort(unique(samp$study_name))) {
  i <- which(samp$study_name == s); cnt <- m[, i, drop = FALSE]; cnt <- cnt[rowSums(cnt) > 0, , drop = FALSE]
  tax <- t(vapply(strsplit(rownames(cnt), "|", fixed = TRUE), function(p) c(sub("^[a-z]__", "", p), rep(NA, 7))[1:7], character(7)))
  tax <- as.data.frame(tax, stringsAsFactors = FALSE); names(tax) <- ranks
  sm <- samp[i, colSums(!is.na(samp[i, , drop = FALSE])) > 0, drop = FALSE]
  out[[s]] <- list(counts = cnt, taxonomy = tax, samples = sm, platform = "shotgun", country = paste(unique(sm$country), collapse = ";"))
}
saveRDS(out, file.path(Sys.getenv("S"), "crc", "shotgun_plain.rds"))
summ <- do.call(rbind, lapply(names(out), function(s) { o <- out[[s]]; d <- colSums(o$counts)
  data.frame(study = s, country = o$country, CRC = sum(o$samples$study_condition == "CRC"), control = sum(o$samples$study_condition == "control"),
             species = nrow(o$counts), genera = length(unique(na.omit(o$taxonomy$genus))), median_reads = round(median(d)),
             age = sum(!is.na(o$samples$age)), sex = sum(!is.na(o$samples$gender)), BMI = sum(!is.na(o$samples$BMI)),
             subjects = length(unique(o$samples$subject_id))) }))
print(summ, row.names = FALSE)
write.csv(summ, file.path(Sys.getenv("S"), "crc", "shotgun_summary.csv"), row.names = FALSE)
