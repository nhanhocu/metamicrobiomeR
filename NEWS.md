# metamicrobiomeR 2.0.0

**Results that differ from version 1.2.** Linear mixed-model p-values in `taxa.compare()` and `pathway.compare()` were
`1.96 * pnorm(-|z|)` in version 1.2, about 2% too small, and are now `2 * pnorm(-|z|)`. Rerunning a 1.x script therefore
gives slightly larger p-values than the 2019 paper (for example 0.0089 in place of 0.0088), with no change in its
conclusions. `meta.taxa()` also pools fewer taxa when few studies have an estimate. Details are below and in
`vignette("migrating-from-1x")`.

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
* `taxa.compare()`, `pathway.compare()`, `taxcomtab.show()`, `metatab.show()`: p-value columns are now found by the literal
  text `Pr(>|t|)`. Before, any column containing "Pr" was taken as a p-value column, so when a covariate name or level
  contained "Pr" (for example a level "Pre" or a variable "Probiotic"), the estimates and standard errors were replaced
  by adjusted values and extra `pval.adjust` columns appeared (GitHub issue #18). Results without "Pr" in a name are unchanged.
* `taxa.compare()`, `pathway.compare()`: taxa and pathways whose names are not syntactic R names (with spaces, semicolons or
  dashes) are now fitted. Before, they failed with an error and were left out (#5).
* `taxcomtab.show()`, `metatab.show()`, `taxa.mean.plot()`: a table holding a single taxonomic level (for example phyla only)
  gave an empty result. It is now shown (#22).
* `taxa.mean.plot(tax.lev = "l2")` keeps phylum-level taxa only, so kingdom-level rows such as `k__archaea` are no longer
  stacked with the phyla (#3).

Changes that keep legacy results the same:

* `taxa.compare()`, `pathway.compare()`: `adjustvar` defaults to `NULL`, so a model without adjustment variables works (#17).
* `taxcomtab.show()`, `taxa.mean.plot()`: `tax.select` accepts several taxa, and `taxa.meansdn()` several variables in
  `othervar`. Before, these stopped with "the condition has length > 1" in R 4.2 and later (#2, #7). `taxa.meansdn()` also
  leaves `sumvar` and `groupvar` out of the abundance variables when `othervar` is given.
* `metatab.show()` uses the study labels when the table has no `pop` column (#19), and warns when the table has no taxa at the
  chosen `tax.lev`, instead of returning an empty table silently (#22).
* `taxa.compare()`, `taxa.filter()`, `taxa.meansdn()`: when no column name contains `k__`, they stop with a clear message,
  which explains how to rename QIIME 2 names starting with `d__` (#14).

* `meta.taxa()` sets `method.tau = "DL"`. The `meta` package now defaults to REML, which changed results and stopped with an
  error on the full `tabsex4` table. The stored 2019 pooled result for `genderMale` is reproduced for all 185 taxa.
* `taxa.compare()`, `pathway.compare()`: the new `n.cyc` argument keeps its default of 20. Non-convergence now gives one warning
  with a count, and `attr(result, "fit_info")` holds the details for each taxon.
* `taxa.compare()`, `pathway.compare()`: a new warning flags taxa with no more non-zero samples than mean-model coefficients.
  Such a taxon passes the prevalence filter in a small study, and its estimate and standard error are not reliable.
* `taxa.compare()`, `pathway.compare()`: a new warning says that GAMLSS standard errors with a subject random effect can be too small.
* `meta.taxa()`: a new warning flags results where one study holds over 90% of the random-effects weight, and
  `attr(result, "dominance")` lists them.
* `meta.taxa()`: a new warning flags estimates with a standard error below 2% of the median of their study, and
  `attr(result, "small_se")` lists them. On the stored 2019 example it flags 23 estimates, all from the USA(NC) cohort,
  including Coprococcus. The pooled results are unchanged. The help page also showed only part of its value section,
  because a percent sign was not escaped.
* Adjusted p-value columns use the full covariate name, such as `pval.adjust.age.sample`. The old truncated name,
  `pval.adjust.sample`, stays as an alias so that scripts do not break.

`vignette("notes-earlier-results")` lists what we found when re-running the 1.x examples.

Other changes:

* Fixed the CRAN note for the global variable `id` in `metatab.show()`, which came from the dplyr upgrade. The package now imports
  only what it uses from dplyr, and no longer attaches `gamlss` when loaded. Scripts that relied on that need `library(gamlss)`.
* `class(x) == "try-error"` is replaced by `inherits()`. `taxtab6` is a plain `data.frame`.
* `microbiomeage()` gains a `seed` argument. Its example and the `read.multi()` example no longer download anything.
* `gtab.3stud` is no longer in the package, because it holds sample-level data obtained from investigators.
* `lmerTest` and other unused dependencies are removed. `gdata` and `matrixStats` are replaced by base R. Old scripts that
  call `gdata::drop.levels()` themselves need `gdata` installed, or base `droplevels()` (see `vignette("migrating-from-1x")`).
* `meta.niceplot()` draws its lines with `linewidth` and `geom_errorbar(orientation = "y")`, so current ggplot2 versions no
  longer print deprecation notes. The plots are unchanged. ggplot2 3.4.0 or later is required.
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
* New engine `bezi_glmmtmb`: the zero-inflated beta model of `bezi_gamlss` on relative abundance, fitted with `glmmTMB`, with an
  optional subject random intercept and prevalence part. It does not need counts. In null simulations, GAMLSS standard errors
  with `random()` were too small for sparse features and for covariates that vary between subjects, so `bezi_gamlss` with a
  subject effect now points to `bezi_glmmtmb` once per session. `fit_study()` now uses `bezi_glmmtmb` by default when the
  model has a subject random intercept and glmmTMB is installed, and `bezi_gamlss` otherwise. Set `engine = "bezi_gamlss"`
  to reproduce version 1.x analyses.
* `fit_study()`: the `min_prevalence` and `min_mean` filters now count only the samples with complete model variables,
  the same samples the model uses. A `bezi_glmmtmb` estimate without a finite standard error is
  marked `converged = FALSE`.
* `fit_study()` no longer fits features with too few non-zero samples: a feature needs more non-zero samples than the mean model
  has coefficients. A warning gives the count, `attr(, "low_support")` lists them, and `min_nonzero` changes the threshold
  (`min_nonzero = 0` fits them all). `attr(, "fit_info")` gains `n_nonzero`.
* Diversity: `alpha_diversity()`, `fit_alpha()`, `fit_beta()` (PERMANOVA) and `pool_beta()`.

## Pooling and meta-analysis

* `pool_effects()`: REML with the ad hoc Hartung-Knapp test by default, heterogeneity, prediction intervals, `max_weight`
  and `top_study`, false discovery rate within rank, part and term. REML, DerSimonian-Laird and Paule-Mandel need no
  extra package; `metafor` is used for EB, ML, HE and SJ. Estimates from different engines or scales are refused.
* REML pooling: when the estimating equation for the between-study variance has several roots, `pool_effects()` now
  picks the one with the highest restricted likelihood. Before, it took the first root found, which in about 1 in 1,000
  random problems was a lower maximum and could change the p-value substantially. Paule-Mandel solves its equation
  exactly, so it can differ slightly from metafor, whose iterations stop earlier. The help page now says which
  estimators use the internal engine and which use metafor.
* When features are pooled from fewer than 5 studies with a Hartung-Knapp test, `pool_effects()` says so once per session
  and suggests `test = "z"` or `pool_sensitivity()` as a sensitivity analysis. The default test is unchanged: in
  simulations where the studies truly differed, the z test gave 9% to 20% false positives and the default 1% to 7%. Version 1.x
  (`meta.taxa()`) pooled with DerSimonian-Laird and the normal test, which is less conservative with few studies.
* `leave_one_out()`, `meta_regress()` with `study_table()`, `pool_sensitivity()` (the same effects under several settings),
  `combine_pvalues()` (directional Stouffer or Fisher with a sign test; works across engines and tools),
  `adjust_hierarchical()` (error control across taxonomic ranks).

## Checks and interpretation

* `check_poolable()` (mixed measured and predicted functions, databases, regions, small studies, taxon overlap),
  `feature_coverage()`, `check_effects_quality()` and `filter_effects()` (implausibly small standard errors,
  non-converged fits, large estimates). `check_effects_quality()` also compares each standard error with the median
  of its own study (`within_ratio = 0.02`), which works with a single study and catches estimates from features with too
  little information. `fit_study()` warns about such estimates, and `pool_effects()` counts them in the new column
  `n_flagged_se` and warns when a pooled result includes one.
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
* `fit_study(engine = "zibb_glmmtmb")` warns when the median library size is above one million, as in shotgun profiles
  or proportions scaled to read numbers. In 11 public shotgun studies the count model converged for only 40% to 72% of
  genera at that depth, against 86% to 100% for 16S read counts. `meta_analyse(engine = "auto")` then uses `lm_asin`.
* `vignette("migrating-from-1x")` explains that the 1.x model reports only the abundance part, and how to report the
  prevalence part as well.
* `meta_analyse(engine = "auto")` uses `lm_clr` with repeated samples (a study with a subject column, unless
  `fit_args = list(longitudinal = FALSE)`). In the benchmark it kept false positives between 4% and 6% in every setting with
  repeated samples, while `zibb_glmmtmb` and `bezi_glmmtmb` exceeded 7.5% with sparse taxa and `zibb_glmmtmb` failed to
  converge in up to a fifth of sparse data sets. Without repeated samples the choice is unchanged.
* `plot(res, type)` draws `"forest"`, `"heatmap"`, `"volcano"`, `"alpha"`, `"beta"` or `"influence"`; `plot_forest()`,
  `plot_heatmap()` and `plot_volcano()` take effect tables from the step-by-step workflow.
* `write_report()` writes a Markdown report with figures; `manifest()` and `write_manifest()` record data and settings.
* `simulate_study()` and `simulate_studies()` for tests and benchmarks.

## Messages and dependencies

* `citation("metamicrobiomeR")` gives the 2019 paper and the package. The README has the GitHub install line and the tags of
  older versions.

* Informational messages say what the package decided or left out (samples with missing values, features below the
  filters, features not pooled, unassigned reads, a subject effect added by default) and stay quiet otherwise;
  `options(metamicrobiomeR.quiet = TRUE)` silences them. Errors name the IDs and values involved and say what to do.
  Third-party console output is captured and reported once. See `dev/MESSAGES.md`.
* The dependencies of version 1.x (`dplyr`, `gamlss`, `gridExtra`, `lme4`, `meta`, `plyr`, `RColorBrewer`, `tidyr`) stay
  required, so the legacy functions work as before. `glmmTMB`, `vegan`, `metafor`, `phyloseq`, `SummarizedExperiment` and
  `biomformat` are optional, and a function that needs one says which and how to install it.
* The help pages of the six legacy datasets cite their sources by DOI (Subramanian et al. 2014 and Ho et al. 2019) in place
  of a link to a laboratory web page that often timed out.
