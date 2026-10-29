# metamicrobiomeR 2.0: specification (v1.0 of the draft)

Version 2.0 is a clean break from 1.x. Two contracts hold the package together: the
**study object** (input) and the **effect table** (what models return and pooling consumes).

## 1. Study object (`mm_study`)

One cohort or dataset processed with one pipeline. Orientation is **features x samples**.

| Component | Rules |
|---|---|
| `assays$counts` | non-negative integers; optional if `relabund` is given |
| `assays$relabund` | proportions in [0, 1]; percent input is converted; optional if `counts` is given |
| `features` | row names are unique feature IDs; optional rank columns kingdom to species; optional `rank` (the level at which a feature is defined; never inferred from the taxonomy) |
| `samples` | one row per sample, row names are the sample IDs |
| `study` | required `id`, `feature_type` (`otu`, `asv`, `taxon`, `function`); recommended `region`, `platform`, `database`, `db_version`, `pipeline`, `measured`, `subject` |

`check_study()` enforces: unique non-missing IDs identical across assays; no missing values; integer
counts; proportions with sums at most 1 within each `rank` group; the subject column exists. Placeholders
for unassigned lineages (`other`, empty tokens, `Unassigned`) are `NA` in the rank columns.

Created by `as_study()` and the importers. `agglomerate()` sums features to a rank and keeps the
unassigned reads in one feature so that totals and proportions are unchanged.

## 2. Effect table (`mm_effects`)

One row per study x feature x model part x term: `study`, `feature`, `rank`, `part`
(`abundance`, `prevalence`, `dispersion`, `overall`), `term`, `estimate`, `se`, `df`, `p_value`, `n`,
`engine`, `effect_scale`, `converged`. `bind_effects()` refuses duplicates and refuses to combine
studies with different effect scales or different engines for the same part and term.

| Engine | Input | Scale |
|---|---|---|
| `zibb_glmmtmb` | counts | `logit_mean` (abundance); `logit_zero_prob` (prevalence, opt-in) |
| `bezi_gamlss` | proportions | `logit_mean` |
| `lm_asin` | proportions | `asin_sqrt_mean` |
| `lm_clr` | counts or proportions | `clr_mean` |
| `fit_alpha()` | any | `sd_units` or `alpha_raw` |

## 3. Pooled table (`mm_pooled`)

`pool_effects()`: random-effects meta-analysis per feature, part and term (`metafor`); columns `k`,
`estimate`, `se`, confidence limits, `statistic`, `p_value`, `tau2`, `I2`, `Q`, `Q_p`, prediction interval,
`max_weight` and `top_study` (influence), `p_adjust` (FDR within the chosen family, default within
rank, part and term). `test = "adhoc"` (default), `"knha"` or `"z"`. `leave_one_out()` and
`meta_regress()` build on it; `adjust_hierarchical()` adds top-down FDR across ranks.

## 4. Not covered

Distance-based results (PERMANOVA) have no estimate and standard error, so they are combined by
p-value (`pool_beta()`), not pooled as effects. Taxonomy harmonisation across databases is by name
(`harmonise_names()`, `agglomerate(key = "name")`) and is lossy for GTDB. No engine corrects for the
compositionality of the data beyond the CLR sensitivity analysis (see `dev/BENCHMARK-PILOT.md`).
