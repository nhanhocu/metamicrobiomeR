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
