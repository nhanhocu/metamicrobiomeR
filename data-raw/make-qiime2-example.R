# Builds inst/extdata/qiime2-example/: SIMULATED data laid out as QIIME 2 writes it, for trying read_qiime2().
# It is made from CohortA of `example_studies`; nothing here is real data.
#   feature-table.tsv    what `biom convert --to-tsv` writes (comment line, "#OTU ID" header, hash-like feature IDs, counts as 5.0)
#   taxonomy.tsv         what `qiime tools export` writes for a taxonomy artifact (Feature ID, Taxon, Confidence; "; " separators)
#   sample-metadata.tsv  a QIIME 2 metadata file with the #q2:types row
# As in real projects: two samples of the table have no metadata, two metadata rows are controls that are not in the table,
# some features are unclassified or classified only to the phylum, and age is declared numeric.
suppressMessages(library(metamicrobiomeR))
set.seed(7)
out <- "inst/extdata/qiime2-example"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
s <- example_studies$CohortA; cn <- mm_counts(s); f <- mm_features(s); smp <- mm_samples(s)
ids <- vapply(seq_len(nrow(cn)), function(i) paste(sprintf("%08x", sample.int(2^30, 4)), collapse = ""), ""); stopifnot(!anyDuplicated(ids))
rownames(cn) <- ids; colnames(cn) <- sub("CohortA_", "", colnames(cn))
writeLines(c("# Constructed from biom file", paste0("#OTU ID\t", paste(colnames(cn), collapse = "\t")),
             apply(cbind(ids, apply(cn, 2, function(x) sprintf("%.1f", x))), 1, paste, collapse = "\t")), file.path(out, "feature-table.tsv"))
lin <- paste0("d__Bacteria; p__", f$phylum, "; g__", f$genus); lin[c(5, 17)] <- "Unassigned"; lin[9] <- paste0("d__Bacteria; p__", f$phylum[9])
write.table(data.frame(`Feature ID` = ids, Taxon = lin, Confidence = round(runif(nrow(cn), 0.7, 1), 3), check.names = FALSE),
            file.path(out, "taxonomy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
keep <- colnames(cn)[-c(1, 2)]
md <- data.frame(`sample-id` = c(keep, "blank1", "blank2"), group = c(as.character(smp$group[-c(1, 2)]), "control", "control"),
                 age = c(smp$age[-c(1, 2)], NA, NA), sex = c(as.character(smp$sex[-c(1, 2)]), NA, NA), check.names = FALSE)
writeLines(c(paste(colnames(md), collapse = "\t"), "#q2:types\tcategorical\tnumeric\tcategorical",
             apply(md, 1, function(r) paste(ifelse(is.na(r), "", r), collapse = "\t"))), file.path(out, "sample-metadata.tsv"))
writeLines(c("Simulated data in the layout QIIME 2 writes. Not real data.", "",
             "feature-table.tsv   the output of `biom convert --to-tsv` (50 samples)",
             "taxonomy.tsv        the output of `qiime tools export` for a taxonomy artifact",
             "sample-metadata.tsv a QIIME 2 metadata file with the #q2:types row (48 samples with metadata; 2 are controls not in the table)",
             "", "Built by data-raw/make-qiime2-example.R. Try:  ?read_qiime2"), file.path(out, "README.txt"))
cat("written:", paste(list.files(out), collapse = ", "), "\n")
