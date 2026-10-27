# metamicrobiomeR 2.0: format specification (draft 0.1)

Status: draft for discussion. Version 2.0 is a clean break: there is no support
for the 1.x wide-table layout. Two contracts hold the package together:
the **study object** (input) and the **effect table** (what models return and
what pooling consumes).

## 1. Study object (`mm_study`)

A study is one cohort or dataset, processed with one pipeline. Orientation is
**features x samples** (as in phyloseq and SummarizedExperiment).

| Component | Type | Rules |
|---|---|---|
| `assays$counts` | numeric matrix, features x samples | Non-negative integers. Optional if `relabund` is given. |
| `assays$relabund` | numeric matrix, same dimnames | Proportions in [0, 1]. Optional if `counts` is given. Percent input is converted on entry. |
| `features` | data.frame, one row per feature | Row names are unique feature IDs (ASV/OTU ID, taxon string, pathway ID). Optional rank columns: kingdom, phylum, class, order, family, genus, species. Optional `rank`: the level at which the feature is defined (a taxonomic rank for agglomerated tables such as QIIME1 L2 to L7; `leaf` for ASV/OTU features). It is never inferred from the taxonomy: a genus-level feature with an unassigned genus is still genus-level. |
| `samples` | data.frame, one row per sample | Row names are sample IDs matching the assay columns. |
| `study` | named list | Required: `id`, `feature_type`. Recommended: `region`, `platform`, `database`, `db_version`, `pipeline`, `measured` (functions only), `subject` (name of the repeated-measures ID column in `samples`). |

`feature_type` is one of `otu`, `asv`, `taxon`, `function`.

Rules enforced by `check_study()`:

* Feature IDs and sample IDs are unique, non-missing, and identical across assays.
* No missing values in assays; counts are non-negative integers.
* Relative abundances are in [0, 1]; sums within each `rank` group may not exceed 1 (all features form one group if `rank` is absent)
  (sums below 1 are allowed: filtering and unassigned reads are normal).
* `samples` has one row per assay column; `study$subject` exists in `samples` if given.
* Warnings (not errors) for samples with zero depth.

Taxonomy lives in rank columns, never in feature names. Placeholders for
unassigned or padded lineages (`other`, empty tokens such as `g__`, `Unassigned`)
are stored as `NA` in the rank columns; importers are responsible for this. Ranks are analysed
separately (rank-aware agglomeration comes in a later release); features at
different ranks are nested, not independent.

## 2. Effect table (`mm_effects`)

One row per study x feature x model part x term. This is the only thing
the pooling layer needs, so any engine (or an adapter for another tool's
output) can feed it.

| Column | Type | Meaning |
|---|---|---|
| `study` | character | Study ID |
| `feature` | character | Feature ID |
| `rank` | character | Feature rank (NA if not applicable) |
| `part` | character | `abundance`, `prevalence`, `dispersion` or `overall` |
| `term` | character | Model term, e.g. `sexMale` |
| `estimate` | numeric | Coefficient |
| `se` | numeric | Standard error (>= 0 or NA) |
| `df` | numeric | Degrees of freedom, NA if not defined |
| `p_value` | numeric | Within-study p-value, NA if not available |
| `n` | integer | Samples used |
| `engine` | character | Engine that produced the estimate |
| `effect_scale` | character | Scale of `estimate`, e.g. `logit_mean`, `log_odds`, `smd` |
| `converged` | logical | Did the fit converge |

Pooling is only valid when every study reports the same `effect_scale` for a
term; `bind_effects()` refuses to mix scales unless told to.

## 3. Out of scope for the first release of the format

Rank agglomeration, import helpers (phyloseq, TreeSummarizedExperiment, BIOM,
DADA2, MetaPhlAn), harmonisation across taxonomies and databases, alpha and beta
diversity outcomes, and function-level metadata beyond `measured`.
