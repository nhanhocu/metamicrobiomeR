# Pilot benchmark of the per-study engines

Status: a small pilot to check the machinery and give early signals. It is **not** a
paper-grade benchmark. Scripts: `data-raw/benchmark-pilot.R`, `data-raw/benchmark-summary.R`.

## Design

One study of 60 independent samples and 12 features, simulated with `simulate_study()`
(Dirichlet-multinomial, depth varying with sd of the log depth 0.8, overdispersion 20, no
structural zeros), 60 replicates per scenario and engine. Every engine tests the group term
for every feature at p < 0.05; the false discovery proportion (FDP) uses Benjamini-Hochberg at 5%
within a replicate. A feature is "differential" if its generating log-abundance parameter has a
group effect. Because the data are compositional, a change in some features also shifts the
*proportions* of the others, which this truth definition counts as false positives.

| Scenario | Differential features | Effect (log scale) |
|---|---|---|
| null | 0 | 0 |
| alt_balanced | 4 (2 up, 2 down) | 1.0 |
| alt_unbalanced | 4, all up | 1.5 |

## Results (60 replicates, 716 or 717 tests per cell)

| Scenario | Engine | False-positive rate (non-differential, p<0.05) | Power | Mean FDP at BH 5% | Converged | s per replicate |
|---|---|---|---|---|---|---|
| null | bezi_gamlss | 0.041 | | | 1.000 | 0.96 |
| null | lm_asin | 0.045 | | | 1.000 | 0.04 |
| null | lm_clr | 0.035 | | | 1.000 | 0.04 |
| null | zibb_glmmtmb | 0.057 | | | 1.000 | 2.61 |
| alt_balanced | bezi_gamlss | 0.264 | 0.703 | 0.312 | 1.000 | 0.94 |
| alt_balanced | lm_asin | 0.291 | 0.874 | 0.306 | 1.000 | 0.04 |
| alt_balanced | lm_clr | 0.101 | 0.916 | 0.098 | 1.000 | 0.04 |
| alt_balanced | zibb_glmmtmb | 0.340 | 0.921 | 0.341 | 0.999 | 2.63 |
| alt_unbalanced | bezi_gamlss | 0.512 | 0.738 | 0.548 | 1.000 | 1.01 |
| alt_unbalanced | lm_asin | 0.606 | 0.850 | 0.548 | 1.000 | 0.04 |
| alt_unbalanced | lm_clr | 0.543 | 0.958 | 0.490 | 1.000 | 0.04 |
| alt_unbalanced | zibb_glmmtmb | 0.696 | 0.879 | 0.572 | 1.000 | 2.31 |

Monte Carlo interval for a 5% rate with about 716 tests: roughly +/- 0.016 (wider in practice,
because tests within a replicate are correlated).

## What this suggests

1. **No type I problem when nothing differs.** All engines are near 5% under the null.
2. **Compositional bias dominates once something differs.** Every engine on proportions or counts
   reports many non-differential features as different, and the realised FDR is far above 5%.
   Only the CLR engine in the balanced scenario comes close (about 10%); with one-directional
   changes it also fails. None of the engines corrects for compositionality, so a reference-based
   or ratio approach is a gap.
3. **Power:** CLR and the count-based model are highest; GAMLSS-BEZI is lowest. The data are not
   BEZI-generated, which may favour the count model.
4. **Convergence at 100 iterations** was essentially complete for every engine.
5. **Cost:** the count model is about 60 times slower than the linear models and 2.5 times
   slower than GAMLSS (single core).

## Limits

12 features, one sample size, one overdispersion, one simulation family, single studies only,
no structural zeros, no repeated measures, no depth difference between groups (the scenario that
should separate count-based from proportion-based engines), and the pooling step is not
benchmarked here. Treat every number above as preliminary.
