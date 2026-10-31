# Cross-platform benchmark: two-stage meta-analysis versus a combined analysis

Script: `data-raw/benchmark-crossplatform.R`; raw results: `dev/benchmark-crossplatform-results.csv`. Preliminary: one scenario.

## Design
Nine simulated replicates. Each has three studies of 50 samples on different platforms: two 16S studies (45 and 50
genera) and a shotgun study (45 genera, 25 times the depth, GTDB-style names such as `Genus20_A`, bracketed names such as
`[Genus05]`), a platform-specific distortion of every genus, and person-to-person variation. Five genera present in all
studies have a true difference (+1.5 or -1.3 on the log-odds scale); one further genus (two studies only) has +1.2.
Compared, with the count model (`zibb_glmmtmb`) and the CLR linear model (`lm_clr`):

* **two-stage**: `meta_analyse()` with the default test (REML, ad hoc Hartung-Knapp) and with the normal test;
* **combined**: `merge_studies()` then one model with the study as a covariate.

## Results (mean over 9 replicates)

| Design | Genera analysed | Ranks true genera on top (effect AUC) | Direction correct | True genera found at FDR 5% (of 5) |
|---|---|---|---|---|
| Two-stage, default test | 50 | 1.00 | 100% | **0** (in all 9 replicates) |
| Two-stage, normal (z) test | 50 | 1.00 | 100% | 5 |
| Combined analysis | 30 | 1.00 | 100% | 5 |

Time per replicate (count model): two-stage about 47 s, combined about 15 s; CLR: 0.4 s and 0.1 s. One core.

## What this supports
* Both approaches recover the true genera by effect size, in the right direction, on every replicate.
* The two-stage analysis uses 50 genera, the combined analysis only the 30 present in every study.
* The two-stage analysis works without any shared table.

## What it does not support
* **It does not show that meta-analysis is more powerful.** With three studies the default two-stage test found none of
  the five true genera in any replicate (its t distribution has two degrees of freedom), while the combined analysis found all five.
  The normal test finds them but ignores the uncertainty of the between-study variance, and is known to be anti-conservative with few studies.
* False positives were not measured. Relative abundances are compositional, so genera without a true effect shift too and a
  "true negative" is ambiguous. Power alone therefore flatters the permissive options.
* One scenario (three studies, 50 samples, large effects, simple platform distortions). No real copy-number bias,
  taxonomy-resolution differences, or repeated measures. More studies would help the two-stage default test.
