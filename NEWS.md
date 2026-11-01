# metamicrobiomeR 2.0.0.9000

Version 2.0 is a redesign built around two data contracts, the **study object** and the **effect table**
(specified in `dev/SPEC.md`).

## Studies and data

* `as_study()`, `check_study()`, the accessors (`mm_counts()`, `mm_relabund()`, `mm_features()`, `mm_samples()`,
  `mm_depth()`, `mm_study_info()`), `study_table()`, `subset_samples()`, `subset_features()`, `agglomerate()` (taxa
  matched by name across databases) and `merge_studies()` (one-stage analysis with the study as a covariate).
* Importers: `read_dada2()`, `read_qiime2()`, `read_metaphlan()`, `read_biom_study()`, `study_from_phyloseq()`,
  `study_from_se()`, `read_function_table()`; exporters `as_phyloseq()` and `as_summarized_experiment()`;
  `effects_from()` brings in results of MaAsLin2, ANCOM-BC2, LinDA or any tool.
* `read_qiime2()` reads QIIME 2 artifacts (`.qza`, for the taxonomy and the feature table) as well as exported text files, and
  honours the `#q2:types` row of the metadata (numeric columns are numbers, categorical columns are factors). All file importers keep the
  samples present in both the table and the sample data and say which were left out. A simulated example in QIIME 2's layout is in
  `inst/extdata/qiime2-example`. `study_from_phyloseq()`, `study_from_se()` and `read_dada2()` leave out taxonomic ranks with no assignment at all,
  so a study from `qiime2R` and one from `read_qiime2()` have the same ranks.
* Taxonomy: `parse_taxonomy()`, `harmonise_names()`, `feature_parents()`. A `Domain` column is read as the kingdom.
* `example_studies`: a simulated example (four case-control studies on three platforms, real genus names, different
  naming conventions, six genera with a known difference).

## Per-study models

* `fit_study()` with the engines `zibb_glmmtmb` (zero-inflated beta-binomial on counts, uses the sequencing depth,
  optional prevalence part), `bezi_gamlss` (the GAMLSS-BEZI model of Ho et al. 2019; default `n_cycles = 100`),
  `lm_asin` and `lm_clr` (compositional sensitivity). Estimates of different engines cannot be combined. Unassigned
  lineages are excluded by default. `cores` fits features in parallel. The linear engines fit all features of a study in one step.
* Diversity: `alpha_diversity()`, `fit_alpha()`, `fit_beta()` (PERMANOVA) and `pool_beta()`.

## Pooling and meta-analysis

* `pool_effects()`: REML with the ad hoc Hartung-Knapp test by default, heterogeneity, prediction intervals, `max_weight`
  and `top_study`, false discovery rate within rank, part and term. REML, DerSimonian-Laird and Paule-Mandel need no
  extra package; `metafor` is used for EB, ML, HE and SJ. Estimates from different engines or scales are refused.
* `leave_one_out()`, `meta_regress()` with `study_table()`, `pool_sensitivity()` (the same effects under several settings),
  `combine_pvalues()` (directional Stouffer or Fisher with a sign test; works across engines and tools),
  `adjust_hierarchical()` (error control across taxonomic ranks).

## Checks and interpretation

* `check_poolable()` (mixed measured and predicted functions, databases, regions, small studies, taxon overlap),
  `feature_coverage()`, `check_effects_quality()` and `filter_effects()` (implausibly small standard errors,
  non-converged fits, large estimates).
* `summary()` of a result describes the between-study variation (I2, tau2, heterogeneity test), whether conclusions
  would hold in a new study (prediction intervals), the most heterogeneous features and what can explain the differences.

## One call, graphs and reports

* `meta_analyse()` runs the whole workflow (checks, per-study models, quality check, pooling, optionally alpha and beta
  diversity) and returns an object with `print()`, `summary()`, `plot()` and `as.data.frame()` methods. Messages and warnings
  are labelled with the study they come from.
* `plot(res, type)` draws `"forest"`, `"heatmap"`, `"volcano"`, `"alpha"`, `"beta"` or `"influence"`; `plot_forest()`,
  `plot_heatmap()` and `plot_volcano()` take effect tables from the step-by-step workflow.
* `write_report()` writes a Markdown report with figures; `manifest()` and `write_manifest()` record data and settings.
* `simulate_study()` and `simulate_studies()` for tests and benchmarks.

## Messages and dependencies

* Informational messages say what the package decided or left out (samples with missing values, features below the
  filters, features not pooled, unassigned reads, a subject effect added by default) and stay quiet otherwise;
  `options(metamicrobiomeR.quiet = TRUE)` silences them. Errors name the IDs and values involved and say what to do.
  Third-party console output is captured and reported once. See `dev/MESSAGES.md`.
* Only `ggplot2` is required. `glmmTMB`, `gamlss`, `lme4`, `vegan`, `metafor`, `phyloseq`, `SummarizedExperiment` and
  `biomformat` are optional, and a function that needs one says which and how to install it.
