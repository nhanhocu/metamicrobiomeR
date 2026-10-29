# metamicrobiomeR (development version)

Version 2.0 is a redesign and a **clean break from 1.x**: the 1.x wide-table
input layout and the dot-case functions (`taxa.compare()`, `meta.taxa()` and
so on) will not be supported. Version 1.3 (branch `legacy-1.x`) is the last
release with the old interface. See `dev/SPEC.md` for the new data contracts.

* New study object (`as_study()`, `check_study()`, accessors) with taxonomy in
  rank columns, optional counts and validated study metadata.
* New effect table (`new_effect_table()`, `check_effects()`, `bind_effects()`)
  as the interface between models and pooling.
* `fit_study()` with the 1.x GAMLSS-BEZI engine (`bezi_gamlss`): fits a model to
  every feature of a study and returns an effect table, with convergence reported
  for every fit and an `n_cycles` option. It reproduces `taxa.compare()` from 1.x
  and all GAMLSS-BEZI results in Tables 3 to 5 of Ho et al. (2019); see
  `dev/VALIDATION.md` for the validation and the findings it produced.
* `pool_effects()`: random-effects meta-analysis of an effect table (REML by
  default, Hartung-Knapp with the ad hoc correction, heterogeneity, prediction
  intervals, FDR within rank). It reproduces the 1.x pooled sex-difference result
  with `method = "DL", test = "z"`. 1.x also "pooled" taxa seen in only one study
  for terms missing from some studies; 2.0 requires estimates from at least two.
* `fit_study(engine = "lm_asin")`: the arcsine-square-root linear (mixed) model of 1.x,
  reproducing its results (the 1.x mixed-model p-values were 2% too small, see 1.3).
* `fit_study()` defaults to `n_cycles = 100` and drops features unassigned at
  their own rank (`exclude_unassigned = TRUE`); both can be set to the 1.x behaviour.

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
