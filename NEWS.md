# metamicrobiomeR (development version)

Version 2.0 is a redesign and a **clean break from 1.x**: the 1.x input format and functions
are gone. Version 1.3 (branch `legacy-1.x`) is the last release with the old interface; see
`vignette("migrating-from-1x")`. The specification is in `dev/SPEC.md`; validation against the
2019 paper in `dev/VALIDATION.md`.

## New

* **Study object and effect table** as the two data contracts (`as_study()`, `check_study()`,
  `new_effect_table()`, `bind_effects()`), with taxonomy in rank columns and validated study metadata.
* **Importers:** `read_dada2()`, `read_qiime2()`, `read_metaphlan()`, `read_biom_study()`,
  `study_from_phyloseq()`, `study_from_se()`, `read_function_table()`; `as_phyloseq()` and
  `as_summarized_experiment()`; `subset_samples()`, `subset_features()`.
* **Taxonomy:** `parse_taxonomy()`, `harmonise_names()`, `agglomerate()`, `feature_parents()`.
* **Engines in `fit_study()`:** `zibb_glmmtmb` (zero-inflated beta-binomial on counts, uses the
  sequencing depth, optional prevalence part), `bezi_gamlss` (the 1.x model; default `n_cycles = 100`),
  `lm_asin`, and `lm_clr` (compositional sensitivity). Estimates of different engines cannot be combined.
* **Pooling:** `pool_effects()` (REML and the ad hoc Hartung-Knapp test by default, heterogeneity,
  prediction intervals, `max_weight` and `top_study`), `leave_one_out()`, `meta_regress()` with
  `study_table()`, `adjust_hierarchical()`, `pool_beta()`.
* **Checks before pooling:** `check_poolable()` (mixed measured/predicted functions, databases, regions, small
  studies, taxon overlap), `feature_coverage()`, `check_effects_quality()` and `filter_effects()` (implausibly
  small standard errors, non-converged fits, large estimates), `pool_sensitivity()` (the same effects pooled under
  several settings).
* **More ways to combine studies:** `combine_pvalues()` (directional Stouffer or Fisher with a sign test; works
  across engines and tools), `merge_studies()` (one-stage analysis with study as a covariate), `meta_regress()`.
* **Results from other tools:** `effects_from()` for MaAsLin2, ANCOM-BC2, LinDA and generic tables.
* **Plots:** `plot_heatmap()` (studies by features with a pooled column) and `plot_volcano()`.
* **Diversity:** `alpha_diversity()`, `fit_alpha()`, `fit_beta()`.
* **Output and reporting:** `plot_forest()`, `manifest()`, `write_manifest()`.
* **Simulation:** `simulate_study()`, `simulate_studies()`.

## Easier to use, and faster

* `meta_analyse()` runs the whole workflow in one call (checks, per-study models, quality check, pooling) and
  returns an object with `print()`, `plot()` and `as.data.frame()` methods. Messages and warnings are labelled with
  the study they come from. The step-by-step functions remain.
* Pooling is 14 to 45 times faster (a lean DL/REML/PM core, tested against metafor; metafor remains the fallback
  for other estimators), `leave_one_out()` and `pool_sensitivity()` follow, `meta_regress()` and `combine_pvalues()`
  are faster, `parse_taxonomy()` is vectorised (identical output, tested against the original), and `fit_study()`
  gains `cores` for parallel fitting. See `dev/PERFORMANCE.md`.
* Help pages are grouped into families and the package help page lists where to start.

* `summary()` of a `meta_analyse()` result describes the between-study variation (I2, tau2, heterogeneity test), whether
  conclusions would hold in a new study (prediction intervals), the most heterogeneous features, and what can explain
  the differences, including batch effects within studies. The printed summary hints when heterogeneity is high.

## Messages

* Informational messages say what the package decided or left out (samples with missing values, features
  below the filters, features not pooled, unassigned reads, a subject effect added by default) and stay quiet
  otherwise; `options(metamicrobiomeR.quiet = TRUE)` silences them. Errors name the IDs and values involved and
  say what to do. Third-party console output (lme4 singular fits, glmmTMB, GAMLSS) is captured and reported once.
  `fit_info` gains `first_warning` and `singular`. Effect and pooled tables print a short header.
  See `dev/MESSAGES.md`.
* New guards: pooling or meta-regression of estimates from different engines or scales now fails with a clear
  error (previously only `bind_effects()` checked); a `Domain` taxonomy column is read as the kingdom.

## Changed from 1.x

* 1.x mixed-model p-values were 2% too small, adjusted p-value columns were truncated at the last
  dot, the GMPR transform returned only missing values, and `meta.taxa()` "pooled" taxa seen in only
  one study; none of this carries over.
* Unassigned lineages are excluded from models by default; FDR is adjusted within a rank.
* The package now has few dependencies: `ggplot2`, `lme4`, `metafor`; the engines and importers use
  suggested packages.

# metamicrobiomeR 1.3

This is a maintenance release. It restores the package to a working state on
current versions of R and its dependencies; there are no new features. A
redesigned version 2.0 is planned.

## Maintainer and licence

* New maintainer: Payton Yau. Nhan Ho remains an author.
* The licence is now GPL (>= 2) (previously GPL-2).

## Bug fixes

* `taxa.compare()` and `pathway.compare()`: p-values for linear mixed models
  (`propmed.rel = "lm"`, `longitudinal = "yes"`) were computed as
  `1.96 * pnorm(-|z|)` instead of `2 * pnorm(-|z|)`, which made them about 2%
  too small. Results from these models will change slightly; please re-run
  analyses that used them. Models fitted with GAMLSS and non-longitudinal
  linear models are not affected.
* `taxa.compare()`, `pathway.compare()` and `taxcomtab.show()`: the adjusted
  p-value columns were named by cutting the covariate name at its last dot, so
  a covariate such as `age.sample` gave `pval.adjust.sample`. Columns now keep
  the full covariate name (`pval.adjust.age.sample`). Covariate names without a
  dot are unaffected.
* `taxa.compare(transform = "gmpr")`: GMPR size factors were `NA` for every
  sample, so every model failed. They are now calculated. GMPR was designed for
  count data on a full community table; on relative abundances or very small
  taxa subsets, samples sharing too few taxa still get `NA` (with a warning).
* `taxtab6` is now a plain `data.frame`. It was stored as an old-format dplyr
  `grouped_df` that raised an error on subsetting with current dplyr.
* `taxa.compare()`, `pathway.compare()` and `alpha.compare()` now use
  `inherits()` to test for failed model fits, as `class(x) == "try-error"`
  errors in R >= 4.2 for objects with more than one class.

## Other changes

* Removed the unused `lmerTest` import.
* Removed ten unused packages from `Suggests` (`RCurl`, `httr`, `repmis`,
  `jsonlite`, `magrittr`, `foreign`, `mgcv`, `reshape2`, `gplots`, `tsibble`).
* Added regression tests for the fixes above.
