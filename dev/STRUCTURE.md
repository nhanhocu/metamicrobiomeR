# How the repository is organised

## The R code: numbered by workflow stage

R requires the `R/` folder to be flat, so the groups are expressed in the file names. The first digit is the stage of an analysis, so a directory
listing reads in the order you would use the package; files load in alphabetical order, which does not matter because no file runs code at load time
that depends on another.

| Digit | Stage | What belongs there |
|---|---|---|
| 0x | Package and helpers | package documentation, message helpers |
| 1x | Studies | the study object, operations on studies (subset, agglomerate, merge), taxonomy |
| 2x | Import | files, Bioconductor objects, results of other tools |
| 3x | Effect table | what models return and pooling consumes |
| 4x | Per-study models | abundance engines, alpha and beta diversity |
| 5x | Pooling | pooling, its engine, meta-regression, sensitivity and p-values, hierarchical control |
| 6x | Checks | checks before and after pooling |
| 7x | Plots | all graphs |
| 8x | Workflow and output | `meta_analyse()`, reports, manifests |
| 9x | Simulation and data | simulated studies, the example data |

| File | What is in it | Exported functions | Tests |
|---|---|---|---|
| `00-package.R` | Package documentation and names used in non-standard evaluation | none (internal) |  |
| `01-messages.R` | Message helpers that keep the wording of errors, warnings and notes consistent (see dev/MESSAGES.md) | none (internal) | `test-01-messages.R` |
| `10-study-object.R` | Study object: the input contract (see dev/SPEC.md) | `as_study()`, `check_study()`, `mm_counts()`, `mm_relabund()`, `mm_depth()`, `mm_features()`, `mm_samples()`, `mm_study_info()`, `study_table()` | `test-10-studies.R` |
| `11-study-operations.R` | Operations on studies: subset, bring to a taxonomic rank, merge | `subset_samples()`, `subset_features()`, `agglomerate()`, `merge_studies()` | `test-10-studies.R` |
| `12-taxonomy.R` | Taxonomy: parse the naming conventions of different databases, harmonise names, find parents | `parse_taxonomy()`, `harmonise_names()`, `feature_parents()` | `test-12-taxonomy.R` |
| `20-import-files.R` | Importers that read files: DADA2, QIIME 2, BIOM, MetaPhlAn and function tables | `read_dada2()`, `read_qiime2()`, `read_metaphlan()`, `read_biom_study()`, `read_function_table()` | `test-20-importers.R` |
| `21-import-objects.R` | Importers and exporters for Bioconductor objects (phyloseq, SummarizedExperiment) | `study_from_phyloseq()`, `study_from_se()`, `as_phyloseq()`, `as_summarized_experiment()` | `test-20-importers.R` |
| `22-import-results.R` | Results of other tools as an effect table | `effects_from()` | `test-20-importers.R` |
| `30-effects-table.R` | Effect table: what models return and what pooling consumes (see dev/SPEC.md) | `new_effect_table()`, `check_effects()`, `bind_effects()` | `test-30-effects.R` |
| `40-model-abundance.R` | Per-study differential abundance: the engines of fit_study() | `fit_study()` | `test-40-models.R` |
| `41-model-diversity.R` | Alpha diversity and beta diversity (PERMANOVA) per study, and combining beta diversity | `alpha_diversity()`, `fit_alpha()`, `fit_beta()`, `pool_beta()` | `test-40-models.R` |
| `50-pool-effects.R` | Pooling effect tables across studies, and leaving studies out | `pool_effects()`, `leave_one_out()` | `test-50-pooling.R` |
| `51-pool-engine.R` | The lean random-effects engine behind pooling | none (internal) | `test-51-pool-engine.R` |
| `52-pool-meta-regression.R` | Meta-regression on study characteristics | `meta_regress()` | `test-50-pooling.R` |
| `53-pool-sensitivity.R` | Pooling settings compared, and combining p-values | `pool_sensitivity()`, `combine_pvalues()` | `test-50-pooling.R` |
| `54-pool-hierarchical.R` | Hierarchical false discovery control across taxonomic ranks | `adjust_hierarchical()` | `test-50-pooling.R` |
| `60-checks.R` | Checks before pooling: compatible studies, taxon overlap, doubtful estimates | `feature_coverage()`, `check_poolable()`, `check_effects_quality()`, `filter_effects()` | `test-60-checks.R` |
| `70-plots.R` | Plots: forest, heatmap, volcano, beta diversity, influence | `plot_forest()`, `plot_heatmap()`, `plot_volcano()` | `test-70-plots.R` |
| `80-meta-analyse.R` | The front door: from a list of studies to pooled results in one call | `meta_analyse()` | `test-80-meta-analyse.R` |
| `81-report-manifest.R` | Reports and provenance manifests | `write_report()`, `manifest()`, `write_manifest()` | `test-80-meta-analyse.R, test-90-simulate-manifest.R` |
| `90-simulate-data.R` | Simulated studies and the example data | `simulate_study()`, `simulate_studies()` | `test-90-simulate-manifest.R, test-91-example-data.R` |

Internal helpers (names starting with a dot) sit next to the code that uses them; a helper used by several files lives with the first stage that needs it.

## Everything else

| Path | What it is |
|---|---|
| `tests/testthat/test-NN-*.R` | tests, numbered like the R files they cover (each file's first line says which) |
| `vignettes/` | `getting-started` (the workflow, step by step) and `cross-platform` (two-stage meta-analysis against a combined analysis) |
| `inst/extdata/qiime2-example/` | a simulated project in QIIME 2's layout (feature table, taxonomy, metadata); built by `data-raw/make-qiime2-example.R` |
| `data/example_studies.rda` | the simulated example; built by `data-raw/make-example-studies.R` |
| `data-raw/` | scripts that build the example data and run the benchmarks |
| `dev/SPEC.md` | the specification of the two data contracts (study object and effect table) |
| `dev/MESSAGES.md` | the style guide for errors, warnings and notes, with a catalogue of them |
| `dev/PERFORMANCE.md` | measured speed and dependency footprint |
| `dev/BENCHMARK-*.md`, `dev/POSITIONING.md` | benchmarks and what can and cannot be claimed |
| `dev/STRUCTURE.md` | this file |

## Adding code

1. Put a new function in the file for its stage; start a new file in that group (next free number) if none fits.
2. Document it with roxygen and give it an `@family` tag (the families match the groups above and the help pages' "See also").
3. Add tests to the matching `test-NN-*.R` file, and a message test if it emits messages (see `dev/MESSAGES.md`).
4. Run `roxygen2::roxygenise()`, the tests, and `R CMD check`.
