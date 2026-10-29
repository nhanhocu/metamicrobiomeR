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
