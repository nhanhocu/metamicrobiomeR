# The 1.x example data: format inventory

**Verdict: the examples that shipped with metamicrobiomeR 1.x are in QIIME 1 format.**
Verified from the raw Bangladesh files in `data-raw/legacy-1.1/extdata/QIIME_outputs/Bangladesh/`;
for the other three studies the evidence is the naming convention only (see the last rows).

| Object or file | Format | Evidence |
|---|---|---|
| `Subramanian_et_al_mapping_file*.txt` | QIIME 1 mapping file | required QIIME 1 columns `#SampleID`, `BarcodeSequence`, `LinkerPrimerSequence` |
| `tax_mapping7/..._L2.txt` to `_L7.txt` | QIIME 1 `summarize_taxa.py` levels 2 to 7, merged with the mapping file | level suffixes, `k__`/`p__`/... lineage names |
| `alpha_div_collated/*.txt` | QIIME 1 `collate_alpha.py` layout: rows are rarefaction depth by iteration, samples are columns | header `sequences per sample`, `iteration`; metric names `observed_species`, `PD_whole_tree` (QIIME 2 uses `observed_features`, `faith_pd`) |
| `picrust/metagenome_at_level1..3_edit.txt` | PICRUSt 1 KEGG pathway tables (samples as columns) | `OTU ID` header, KEGG levels 1 to 3 |
| `taxtab6` (CRAN) and `taxtab6.rm[[1..6]]` | wide table, samples as rows; cumulative taxa tables down to phylum, class, order, family, genus, species | names such as `k__bacteria.p__firmicutes.c__clostridia...`, lower-cased and dot-separated, unassigned levels padded with `.other` |
| `gtab.3stud` (Haiti, CA/FL, NC) | same wide layout; only taxa, `x.sampleid` and `age.sample` | same naming convention; the processing of these three studies is not documented in the repository |
| `alphadat` | list of four QIIME 1 collated alpha-diversity tables | as above |
| `kegg.12` | PICRUSt 1 KEGG level 2 and 3, samples as rows | as above |
| stored results (`taxacom*`, `metab.sex`, ...) | outputs of the 1.x functions | wide `Estimate.<term>` columns |

What the paper documents for Bangladesh (Table 6): V4, Illumina MiSeq, open-reference OTU
picking with UCLUST at 97% against Greengenes 13.8. Nothing in the repository says how the
Haiti, CA/FL and NC tables were produced beyond the paper's Table 6 (V4 MiSeq for Haiti and
CA/FL, V1-V2 Roche 454 for NC).

## What is not in the repository

* No counts, only relative abundances (so the count-based model cannot be demonstrated on these data).
* The sample data (gender, feeding, subject IDs) of Haiti, CA/FL and NC: they came directly from
  the study investigators. Only Bangladesh can be refitted from raw data.
* Anything on ASVs, SILVA or GTDB.

## Consequences for 2.0

1. The 1.x examples and the old vignette cannot be carried over: 2.0 does not read this
   format. The conversion helpers in `data-raw/legacy-helpers.R` exist only to validate 2.0
   against the 2019 results.
2. 2.0 needs new-format example data. The Bangladesh data are public and can be converted, but
   they are proportions of 2013-era Greengenes OTUs; a count-based example and an ASV/SILVA/GTDB
   example must come from other public datasets (to be shortlisted).
3. Because every 1.x example is OTU-level and Greengenes-based, taxonomy harmonisation across
   databases cannot be tested with them.
