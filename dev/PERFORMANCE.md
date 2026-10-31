# Performance

Measured on one core, R 4.3.3, 5,000 features x 5 studies of effect estimates (`data-raw/` scripts were ad hoc;
the numbers below are from the profile run when this was written).

| Task | Before | After |
|---|---|---|
| `pool_effects()` (REML + ad hoc Hartung-Knapp) | 37.2 s | 2.6 s |
| `pool_effects()` (DL, normal test) | 28.3 s | 0.6 s |
| `leave_one_out()` (5 poolings) | 162.8 s | 10.4 s |
| `pool_sensitivity()` (4 settings) | about 150 s | 8.8 s |
| `meta_regress()` (1,000 features) | 14.5 s | 2.6 s |
| `meta_regress()` (5,000 features) | about 72 s | 15.6 s |
| `parse_taxonomy()` (20,000 strings) | 4.5 s | see NEWS |

Already fast and left alone: `agglomerate()` (0.2 s for 20,000 ASVs), `merge_studies()`, `feature_parents()`,
`adjust_hierarchical()`, `check_effects_quality()`.

## What changed

* **Fast pooling core** (`R/rma_fast.R`) for the DL, REML and Paule-Mandel estimators, intercept-only and with
  moderators, and the z, Hartung-Knapp and ad hoc tests. metafor stays as the fallback for EB, ML, HE and SJ.
* **No per-feature scans:** the feature loops used `effects[effects$key == key, ]`, a scan of the whole table per
  feature (quadratic). They now split once.
* **`fit_study(cores = )`** fits features in parallel by forking (Linux and macOS). The count model and GAMLSS take
  about a second per feature, so this is where the time goes in real analyses. It was only tested for identical results on a
  one-core machine, so a speed-up has not been measured.

## Equivalence with metafor (tests in `tests/testthat/test-rma-fast.R`)

* DL: agrees with `metafor::rma.uni()` to about 1e-15; REML (run with a tight convergence threshold): about 1e-10;
  estimates, standard errors, intervals, p-values, tau2, I2, Q, prediction intervals, and for meta-regression the
  coefficients and the omnibus test, for z, knha and adhoc tests, with 1 to 3 coefficients, on hundreds of random cases.
* Paule-Mandel: metafor solves its equation less precisely (its root finder keeps a loose tolerance even when a tight
  one is requested), so results differ from metafor by up to 5e-5. The fast core satisfies the Paule-Mandel equation to 4e-10.
* Against metafor's *default* settings the fast core agrees to about four digits; REML and PM results can differ in
  the fourth digit because metafor stops iterating at a threshold of 1e-5.
* The 2019 pooled sex-difference result is reproduced to about 1e-13 with the fast core.

## Lighter and faster, measured later

* **Hard dependencies:** `lme4` moved from Imports to Suggests. It is used in one place (the subject random effect of the linear
  engines). The packages pulled in by the hard dependencies fell from 57 to 33 (recursive count from the installed package
  database, excluding base packages); 24 were reachable only through `lme4`. Verified in a library without `lme4`:
  the package loads, `meta_analyse()` runs with the linear engine, and a subject random effect gives an install hint.
  (`glmmTMB` itself needs `lme4`, so the count model is unavailable there and `engine = "auto"` falls back to the linear model.)
* **Test suite:** 21 files merged into 13, organised by function family, with no test removed (112 blocks before and after, 1,159 expectations
  passing before and after). The whole suite takes about 80 seconds on one core, mostly the count model (`test-models.R`, about 30 s)
  and `meta_analyse()` (about 20 s). The suite was not timed before the merge, so no speed-up is claimed.
* **Progress lines** (`Study 2 of 4: CohortB`) and an engine note (the count model takes about a second per feature and study)
  make a long run less opaque.

## Second pass (measured on one core, same scripts as above)

| Task | Before this pass | After |
|---|---|---|
| `fit_study()`, `lm_asin`, 2,000 features x 100 samples | 5.4-6.3 s | 0.08 s |
| the same, 10,000 features | (about 30 s) | 0.29 s |
| `pool_effects()`, 5,000 features x 5 studies | 2.6 s | 0.95 s |
| `pool_sensitivity()` (4 settings) | 8.8 s | 3.8 s |
| `leave_one_out()` (5 poolings) | 10.4 s | 4.8 s |

* **Why the linear engines were slow:** profiling showed 78% of the time in `data.frame()` calls (two or three small data frames per
  feature), not in the model fits. The batch path fits every feature with one QR factorisation and builds the result tables in one step. Features
  that `lm()` would warn about (a perfect fit), data with missing values, rank-deficient designs and one-level variables are left to the
  per-feature path, so the output is the same row for row (tests compare both paths, `options(metamicrobiomeR.batch_lm = FALSE)` switches the fast path off).
* **Pooling:** half the time was generic matrix algebra for a one-number problem; the intercept-only case now uses plain arithmetic (tested against the matrix version to 1e-9 and against `metafor`).
* **Dependencies:** `metafor` and `lme4` are optional. Footprints (recursive, installed package database): `ggplot2` 28 packages, `metafor` 8,
  both together 33 (they share three), and `lme4` 24 more; the hard dependencies went from 57 packages to 28 (`ggplot2` and its dependencies).
  `ggplot2` was made optional for a while and checked in a library without it, `metafor`, `lme4` and `glmmTMB` (a full analysis ran and loaded no package
  beyond base R), and then made required again because the plots are part of every result.
* **Not done:** the count model (about a second per feature) and GAMLSS are limited by the model fits themselves; `cores` has only been
  tested for identical results on one core. The number of exports (53) is unchanged.
