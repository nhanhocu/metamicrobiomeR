# Positioning

## The claim we can defend
> Meta-analysis for microbiome studies that do not share a platform, a pipeline or a database.

It rests on five things that exist and are tested:
1. Each study is analysed on its own terms; only the **estimates** must be comparable (no shared feature table).
2. Taxon names are matched across databases (`[Eubacterium]`, `Blautia_A`), and `check_poolable()` says what can and cannot be pooled.
3. Honesty tools: implausible-estimate filter, prediction intervals, a comparison of pooling settings, a heterogeneity summary.
4. The combined analysis is available too (`merge_studies()`), so both can be reported and compared.
5. Results from other tools can be pooled (`effects_from()`), so it complements them.

## Evidence we have, and its limits
| Evidence | Limit |
|---|---|
| Fast pooling matches `metafor` (DL, REML to about 1e-10; Paule-Mandel checked against its own equation) | Equivalence, not a new method |
| Cross-platform mock benchmark (`dev/BENCHMARK-CROSSPLATFORM.md`) | One simulated scenario; no false-positive measure |
| Example data set and `vignette("cross-platform")` | Simulated |

## What we should not say yet
* "The best" or "solves platform differences". Platform differences remain; they show up as heterogeneity. The benchmark does
  **not** show higher power than a combined analysis (the opposite, at three studies).
* Anything about other packages' features until checked against their current documentation.
* That alpha and beta diversity are validated: they are not.

## To earn a stronger claim
1. A larger benchmark: more studies, several effect sizes, false positives measured with a compositional-aware definition of truth, repeated measures.
2. One real multi-platform public data set (licence permitting) as a worked example.
3. A feature comparison with other packages, checked against their **current** documentation. To check, at least: MMUPHin (as far as we
   recall, it offers batch adjustment and meta-analysis), `metafor` (general meta-analysis), MaAsLin2/3, ANCOM-BC2, LinDA (differential abundance).
4. A decision on the default pooling test with few studies.

## Wording
* Do: "works without a common platform", "shows how much your studies disagree", "complements a combined analysis".
* Do not: "more powerful", "corrects batch effects", "validated" (for the diversity functions).
