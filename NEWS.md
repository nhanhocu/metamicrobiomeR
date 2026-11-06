# metamicrobiomeR 2.0.0.9000

## Legacy functions (version 1.x)

The 12 functions and 6 of the 7 datasets of 1.x stay available, with the same arguments and output shape, so existing
scripts keep working. They receive bug fixes only. New work should start with `meta_analyse()`, and
`legacy_to_study()` and `vignette("migrating-from-1x")` help you move across.

Changes that can alter legacy results:

* `taxa.compare()`, `pathway.compare()`: p-values for linear mixed models were `1.96 * pnorm(-|z|)` and are now
  `2 * pnorm(-|z|)`, so they rise by about 2%. GAMLSS and non-longitudinal linear models are not affected.
* `meta.taxa()`: a taxon is pooled for a term only if at least two studies, and at least `percent.meta` of the studies, have an
  estimate for that term. It used to count rows, so some taxa were "pooled" from one study or none. The `tabsex4` result for
  `genderMale` is unchanged. For `bfNo_BF`, 86 of 185 taxa are no longer pooled, and for `bfNon_exclusiveBF` and `age.sample`, 1 each.
* `taxa.compare(transform = "logit")`: taxa with zeros are left out with a clear warning, and it stops if every taxon has zeros.
* `taxa.compare(transform = "gmpr")`: size factors are now calculated. They were `NA` for every sample.

Changes that keep legacy results the same:

* `meta.taxa()` sets `method.tau = "DL"`. The `meta` package now defaults to REML, which changed results and stopped with an
  error on the full `tabsex4` table. The stored 2019 pooled result for `genderMale` is reproduced for all 185 taxa.
* `taxa.compare()`, `pathway.compare()`: the new `n.cyc` argument keeps its default of 20. Non-convergence now gives one warning
  with a count, and `attr(result, "fit_info")` holds the details for each taxon.
* `taxa.compare()`, `pathway.compare()`: a new warning flags taxa with no more non-zero samples than mean-model coefficients.
  Such a taxon passes the prevalence filter in a small study, and its estimate and standard error are not reliable.
* `taxa.compare()`, `pathway.compare()`: a new warning says that GAMLSS standard errors with a subject random effect can be too small.
* `meta.taxa()`: a new warning flags results where one study holds over 90% of the random-effects weight, and
  `attr(result, "dominance")` lists them.
* Adjusted p-value columns use the full covariate name, such as `pval.adjust.age.sample`. The old truncated name,
  `pval.adjust.sample`, stays as an alias so that scripts do not break.

`vignette("notes-earlier-results")` lists what we found when re-running the 1.x examples.

Other changes:

* Fixed the CRAN note for the global variable `id` in `metatab.show()`, which came from the dplyr upgrade. The package now imports
  only what it uses from dplyr, and no longer attaches `gamlss` when loaded. Scripts that relied on that need `library(gamlss)`.
* `class(x) == "try-error"` is replaced by `inherits()`. `taxtab6` is a plain `data.frame`.
* `microbiomeage()` gains a `seed` argument. Its example and the `read.multi()` example no longer download anything.
* `gtab.3stud` is no longer in the package, because it holds sample-level data obtained from investigators.
* `lmerTest` and other unused dependencies are removed. `gdata` and `matrixStats` are replaced by base R.
* The licence is GPL-3, because the GMPR code in `taxa.compare()` comes from a package licensed GPL-3.
* `fit_study(engine = "bezi_gamlss")` no longer stops on a sparse feature with aliased coefficients.
* New maintainer: Payton Yau. Nhan Ho remains an author.

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
* `meta_analyse()` drops doubtful estimates before pooling by default (`quality = "filter"`: fits that did not converge, and standard errors far
  smaller than in the other studies) and says how many; `quality = "warn"` keeps them. Such fits often have missing or near-zero standard errors, and one
  of them can otherwise decide the pooled result. Pooling skips features whose estimates are numerically degenerate instead of failing, and says how many.
  When many fits do not converge, the warning says whether it is the rare features that fail, with the numbers, and suggests a higher `min_prevalence`
  (on soil data with 15 to 36 samples per study, the count model converged for 12% of the genera present in 10% of samples or fewer and for 99% of those present in more than 75%).
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
