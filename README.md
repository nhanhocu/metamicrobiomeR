# metamicrobiomeR

Analyse microbiome studies one at a time and **pool the results across studies**.
Originally described in Ho et al. (2019), *BMC Bioinformatics* 20:188,
<https://doi.org/10.1186/s12859-019-2744-2>.

> **Status: version 2.0 is a pre-release and is not on CRAN.** It is a redesign and is **not
> compatible with 1.x**; the API may still change. Version 1.3 (branch `legacy-1.x`) is the last
> release with the old interface.

## Workflow

```r
library(metamicrobiomeR)

# 1. studies: counts or proportions, taxonomy in rank columns, study metadata
study <- read_dada2(seqtab, taxa, samples, list(id = "cohort1", region = "V4", database = "SILVA 138"))
genus <- agglomerate(study, "genus")

# 2. one model per feature -> effect table
eff <- fit_study(genus, ~ group + age, engine = "zibb_glmmtmb")

# 3. pool across studies
bound  <- bind_effects(eff1, eff2, eff3)         # refuses to mix engines or scales
pooled <- pool_effects(bound, term = "groupb")   # REML + ad hoc Hartung-Knapp, FDR within rank

# 4. look, check influence, report
plot_forest(bound, pooled, feature = "g__Blautia", term = "groupb")
leave_one_out(bound, feature = "g__Blautia", term = "groupb")
write_manifest(genus, "cohort1-manifest.txt")
```

See `vignette("getting-started")` for a runnable example and `vignette("migrating-from-1x")`
for the changes from 1.x.

## What is in the package

| Area | Functions |
|---|---|
| Studies | `as_study()`, `check_study()`, `subset_samples()`, `subset_features()`, `study_table()`, `agglomerate()` |
| Import | `read_dada2()`, `read_qiime2()`, `read_metaphlan()`, `read_biom_study()`, `study_from_phyloseq()`, `study_from_se()`, `read_function_table()`; export with `as_phyloseq()`, `as_summarized_experiment()` |
| Taxonomy | `parse_taxonomy()`, `harmonise_names()`, `feature_parents()` |
| Per-study models | `fit_study()` with engines `zibb_glmmtmb` (counts, uses depth, optional prevalence part), `bezi_gamlss` (the 1.x model), `lm_asin`, `lm_clr`; `fit_alpha()`, `fit_beta()` |
| Pooling | `bind_effects()`, `pool_effects()`, `leave_one_out()`, `meta_regress()`, `pool_beta()`, `adjust_hierarchical()` |
| Output | `plot_forest()`, `manifest()`, `write_manifest()` |
| Simulation | `simulate_study()`, `simulate_studies()` |

The two data contracts (the study object and the effect table) are specified in
[`dev/SPEC.md`](dev/SPEC.md).

## Install

```r
remotes::install_github("paytonyau/metamicrobiomeR@v2-dev")     # 2.0 development version
remotes::install_github("paytonyau/metamicrobiomeR@legacy-1.x") # last 1.x release
```

Engines and importers use optional packages (`glmmTMB`, `gamlss`, `vegan`, `phyloseq`,
`SummarizedExperiment`, `biomformat`); a clear error says which one is missing.

## Validation and benchmark

* [`dev/VALIDATION.md`](dev/VALIDATION.md): the 2.0 engines reproduce the numbers of Tables 3 to 5
  of the 2019 paper and its stored results, and what that exercise found (non-converged fits,
  errors in the printed tables, a headline meta-analysis result carried by a single study).
* [`dev/BENCHMARK-PILOT.md`](dev/BENCHMARK-PILOT.md): a small pilot benchmark of the engines on
  simulated data (preliminary).
* [`dev/LEGACY-DATA.md`](dev/LEGACY-DATA.md): the format of the 1.x example data.
* `data-raw/legacy-1.1/` holds the raw Bangladesh QIIME outputs and the 1.x result objects from the
  original repository, used for validation.

## Authors and licence

Nhan Ho (author), Payton Yau (maintainer). GPL (>= 2).
