suppressMessages(library(metamicrobiomeR))
S <- Sys.getenv("S"); a <- readRDS(file.path(S,"crc","amplicon_plain.rds")); s <- readRDS(file.path(S,"crc","shotgun_plain.rds"))
syn <- c("Escherichia/Shigella" = "Escherichia", "Escherichia_Shigella" = "Escherichia", Clostridium_sensu_stricto = "Clostridium", Ruminococcus2 = "Ruminococcus")
clean <- function(x) trimws(gsub("\u00a0", " ", as.character(x)))
num <- function(x) suppressWarnings(as.numeric(clean(x)))
sexf <- function(x) { x <- toupper(substr(clean(x), 1, 1)); factor(ifelse(x %in% c("M","F"), x, NA), levels = c("F","M")) }
mk <- function(cnt, tax, smp, id, platform, region, pipeline, database, country) {
  storage.mode(cnt) <- "double"; cnt <- round(cnt); storage.mode(cnt) <- "integer"
  rownames(tax) <- rownames(cnt) <- make.unique(rownames(cnt))
  rownames(smp) <- colnames(cnt)
  st <- as_study(counts = cnt, taxonomy = tax, samples = smp,
                 study = list(id = id, feature_type = if (platform == "shotgun") "taxon" else "otu", platform = platform, region = region,
                              pipeline = pipeline, database = database, country = country))
  agglomerate(st, "genus", key = "name", synonyms = syn)
}
studies <- list(); info <- list()
for (k in names(s)) {
  o <- s[[k]]; sm <- o$samples; keep <- sm$study_condition %in% c("CRC", "control")
  smp <- data.frame(group = factor(ifelse(sm$study_condition[keep] == "CRC", "CRC", "control"), levels = c("control", "CRC")),
                    age = num(sm$age[keep]), sex = sexf(sm$gender[keep]), row.names = rownames(sm)[keep])
  studies[[k]] <- mk(o$counts[, keep], o$taxonomy, smp, k, "shotgun", "WGS", "MetaPhlAn 3 (curatedMetagenomicData)", "mpa_v30", o$country)
}
lab <- c(crc_baxter = "BaxterNT_2016_16S", crc_zackular = "ZackularJP_2014_16S", crc_zhao = "WangT_2012_16S", crc_xiang = "ChenW_2012_16S", crc_zeller = "ZellerG_2014_16S")
reg <- c(crc_baxter = "V4", crc_zackular = "V4", crc_zhao = "V3", crc_xiang = "V1-V3", crc_zeller = "V4")
ctry <- c(crc_baxter = "USA;CAN", crc_zackular = "USA", crc_zhao = "CHN", crc_xiang = "CHN", crc_zeller = "FRA")
for (k in names(a)) {
  o <- a[[k]]; sm <- o$samples; keep <- sm$DiseaseState %in% c("CRC", "H")
  age <- switch(k, crc_baxter = sm$Age_s, crc_zackular = sm$age, crc_zeller = sm[["Age (years)"]], NA)
  sex <- switch(k, crc_baxter = sm$Gender_s, crc_zackular = sm$sex, crc_zeller = sm$Gender, NA)
  smp <- data.frame(group = factor(ifelse(sm$DiseaseState[keep] == "CRC", "CRC", "control"), levels = c("control", "CRC")),
                    age = if (length(age) == 1) NA_real_ else num(age[keep]), sex = if (length(sex) == 1) factor(NA, levels = c("F","M")) else sexf(sex[keep]),
                    row.names = rownames(sm)[keep])
  if (k == "crc_zeller") smp$sample_id <- sm[["Sample ID"]][keep]
  studies[[lab[[k]]]] <- mk(o$counts[, keep], o$taxonomy, smp, lab[[k]], "16S", reg[[k]], "MicrobiomeHD (USEARCH, 100% OTUs)", "RDP", ctry[[k]])
}
saveRDS(studies, file.path(S, "crc", "crc_studies_genus.rds"))
tab <- do.call(rbind, lapply(studies, function(x) { sm <- mm_samples(x); d <- mm_depth(x)
  data.frame(study = x$study$id, platform = x$study$platform, region = x$study$region, country = x$study$country,
             CRC = sum(sm$group == "CRC"), control = sum(sm$group == "control"), genera = nrow(mm_counts(x)) - 1,
             median_depth = round(median(d)), age_sex = sum(!is.na(sm$age) & !is.na(sm$sex))) }))
print(tab, row.names = FALSE); write.csv(tab, file.path(S, "crc", "crc_studies_table.csv"), row.names = FALSE)
