# One-off validation script (not part of the package API).
# Converts the four 1.x bundled studies into 2.0 study objects to test that the
# new format can represent real data, and that old QIIME1-style names parse
# and round-trip. Study metadata from Table 6 of Ho et al. 2019.
source("data-raw/legacy-helpers.R")

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
