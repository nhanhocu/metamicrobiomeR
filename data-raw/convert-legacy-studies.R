# One-off validation script (not part of the package API).
# Converts the four 1.x bundled studies into 2.0 study objects to test that the
# new format can represent real data, and that old QIIME1-style names parse
# and round-trip. Study metadata from Table 6 of Ho et al. 2019.
library(metamicrobiomeR)

codes <- c(k = "kingdom", p = "phylum", c = "class", o = "order", f = "family", g = "genus", s = "species")

parse_legacy_names <- function(nms) {
  # Legacy names look like k__bacteria.p__firmicutes.c__clostridia.o__x.f__y.other
  # A trailing ".other" (repeated) pads a lineage that is unassigned at deeper
  # levels, so each ".other" adds one level. Padding is stored as NA.
  parts <- strsplit(nms, "\\.(?=[a-z]__)", perl = TRUE)
  tax <- matrix(NA_character_, length(nms), length(codes), dimnames = list(nms, unname(codes)))
  ntok <- extra <- integer(length(nms)); standard <- logical(length(nms))
  for (i in seq_along(parts)) {
    tk <- parts[[i]]
    last <- tk[length(tk)]
    pad <- regmatches(last, regexpr("(\\.other)+$", last))
    n_pad <- if (length(pad)) lengths(regmatches(pad, gregexpr("\\.other", pad))) else 0L
    tk[length(tk)] <- sub("(\\.other)+$", "", last)
    cd <- substr(tk, 1, 1)
    standard[i] <- all(substr(tk, 2, 3) == "__") && identical(unname(cd), names(codes)[seq_along(cd)]) &&
      length(tk) + n_pad <= length(codes)
    if (!standard[i]) next
    ntok[i] <- length(tk); extra[i] <- n_pad
    val <- substring(tk, 4)
    tax[i, seq_along(tk)] <- ifelse(nzchar(val), val, NA_character_)
  }
  list(tax = as.data.frame(tax, stringsAsFactors = FALSE), ntok = ntok, extra = extra,
       level = ntok + extra, standard = standard)
}

rebuild_names <- function(tax, ntok, extra) {
  vapply(seq_len(nrow(tax)), function(i) {
    if (ntok[i] == 0) return(NA_character_)
    v <- unlist(tax[i, seq_len(ntok[i])], use.names = FALSE); v[is.na(v)] <- ""
    paste0(paste0(names(codes)[seq_len(ntok[i])], "__", v, collapse = "."), strrep(".other", extra[i]))
  }, character(1))
}

legacy_to_study <- function(tab, study, id_col = "x.sampleid") {
  tab <- as.data.frame(tab)
  taxcols <- grep("^k__", colnames(tab))
  p <- parse_legacy_names(colnames(tab)[taxcols])
  # 'rank' is the level of the table (tokens plus '.other' padding), not the deepest assigned taxon
  p$tax$rank <- unname(codes)[p$level]
  ids <- as.character(tab[[id_col]])
  dup_ids <- anyDuplicated(ids) > 0
  ids <- make.unique(ids)
  m <- t(as.matrix(tab[, taxcols, drop = FALSE])); colnames(m) <- ids
  samp <- tab[, -taxcols, drop = FALSE]; rownames(samp) <- ids
  s <- as_study(relabund = m, taxonomy = p$tax, samples = samp, study = study)
  list(study = s,
       qc = data.frame(
         id = study$id, samples = ncol(m), features = nrow(m),
         nonstandard_names = sum(!p$standard),
         roundtrip_ok = sum(rebuild_names(p$tax, p$ntok, p$extra) == colnames(tab)[taxcols], na.rm = TRUE),
         duplicated_sample_ids = dup_ids))
}

data(taxtab6); data(gtab.3stud)
res <- list(
  legacy_to_study(taxtab6, list(id = "bangladesh", feature_type = "otu", region = "V4",
                                platform = "Illumina MiSeq", database = "Greengenes", db_version = "13.8",
                                subject = "personid")),
  legacy_to_study(gtab.3stud$haiti, list(id = "haiti", feature_type = "otu", region = "V4", platform = "Illumina MiSeq")),
  legacy_to_study(gtab.3stud$ca_fl, list(id = "usa_ca_fl", feature_type = "otu", region = "V4", platform = "Illumina MiSeq")),
  legacy_to_study(gtab.3stud$nc, list(id = "usa_nc", feature_type = "otu", region = "V1-V2", platform = "Roche 454"))
)
qc <- do.call(rbind, lapply(res, `[[`, "qc"))
print(qc, row.names = FALSE)
for (r in res) {
  f <- mm_features(r$study)
  cat(sprintf("\n%s: features per rank: %s", r$study$study$id,
              paste(names(table(f$rank)), table(f$rank), sep = "=", collapse = ", ")))
}
cat("\n")
invisible(res)
