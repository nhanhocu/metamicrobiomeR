# Validation of the 2.0 engine against Ho et al. (2019)

Status: first gate of the 2.0 plan, extended to the pooling step (section 4) and the LMAS engine (section 5). The 2.0 `fit_study()` engine `bezi_gamlss`
(the 1.x GAMLSS-BEZI model) is compared with (a) the result objects stored in
the original GitHub repository (`data-raw/legacy-1.1/results/`) and (b) the numbers
printed in Tables 3 to 5 of the paper. Scripts: `data-raw/validate-2019-tables.R`, `data-raw/validate-2019-meta.R`, `data-raw/validate-2019-lmas.R`, `data-raw/refit-sex-bangladesh.R`,
`data-raw/convergence-sensitivity.R`, `data-raw/convergence-all-taxa.R`.
Environment: R 4.3.3, gamlss 5.5-5 (the paper used gamlss 5.0-5).

## 1. Reproduction

| Comparison | Result |
|---|---|
| New engine vs numbers printed in Tables 3, 4, 5 (GAMLSS-BEZI columns): estimate, both 95% limits, p, FDR for 16 phylum rows | **80 of 80 match** (to the printed 2 or 4 decimals) |
| New engine vs stored 1.x results, phylum level, 4 examples | max absolute difference about 1e-15 |
| New engine vs stored 1.x results, all 124 taxa of Example 1 | same 124 taxa; max absolute difference 2e-11 |
| New engine vs 1.x `taxa.compare()`, same data (unit test) | equal to 1e-8 |

Not yet reproduced: the simulation studies and the type I error analysis on real data.

## 2. Findings that matter for the paper and for 2.0

**F1. Convergence.** GAMLSS flagged the fit as not converged within its default 20
iterations for 109 of 124 taxa in Example 1 (all ranks, longitudinal model). For the
four main phyla shown in the tables, 14 of 16 fits converged; the two that did not
(Bacteroidetes in Examples 1 and 2) converge with more iterations and move by at most
0.006 in the estimate and 0.005 in the p-value, with no change in significance.
For a random sample of 25 of the 109 non-converged taxa, all converged with 300
iterations, but the results changed noticeably: median |change in p| 0.03 to 0.25
depending on the term, and significance at 0.05 flipped in 3 to 5 of 25 fits per term.
The published phylum-level conclusions are robust; lower-rank results under the 1.x
defaults are not. 2.0 should report convergence for every fit and use a higher
iteration limit by default.

**F2. Example 3 was not fitted as described.** The stored results (and the supplementary
file) use `longitudinal = "no"`, ignoring that the samples come from few infants (377
samples from 30 infants for EBF <= 2 months; 296 from 20 for EBF > 2 months), while the
paper describes longitudinal data with subject random effects. Fitting the longitudinal
model changes the numbers (e.g. Actinobacteria, EBF <= 2 months: estimate -0.73 to -0.60,
p 0.0003 to 0.0012; EBF > 2 months: Actinobacteria p 0.92 to 0.42) but not the
qualitative conclusions (differences for Actinobacteria and Firmicutes only in the
short-EBF group).

**F3. Only the mean part of the model depends on the covariates.** The GAMLSS-BEZI
summary stacks the mean (mu), dispersion (sigma) and zero-inflation (nu) parts; 1.x keeps
the non-intercept rows, and sigma and nu are intercept-only. A group effect on the
probability of absence is therefore never estimated, although the paper describes both
parts as allowing covariates. 2.0 should model the prevalence part explicitly.

**F4. The FDR in the tables is not the FDR the function returns.** The paper's FDR values
are Benjamini-Hochberg adjustments within the four main phyla; the `pval.adjust.*` columns
of `taxa.compare()` adjust across all taxa of all ranks (e.g. Actinobacteria, Example 1:
0.0571 stored vs 0.0166 in the paper). 2.0 should make the family of tests explicit
(per rank, hierarchical).

**F5. Other.** Stored adjusted p-value columns carry the truncated names fixed in 1.3
(`pval.adjust.sample`, `pval.adjust.food5>5 months`). The supplementary tutorial shows
different LMAS numbers for Example 1 (Actinobacteria -0.11, p 0.0066) from the paper's
Table 3 (-0.13, p 0.0088); to be resolved when LMAS is reproduced.

## 3. Decisions taken from these findings (2.0 defaults)

* `fit_study(n_cycles = 100)`, with a warning for any fit that still does not converge.
  In a test, none of 3 fits converged at 20 iterations and all converged at 100, with
  results identical to 300 (earlier: 25 of 25 converged at 300). `n_cycles = 20` reproduces 1.x.
* `fit_study(exclude_unassigned = TRUE)`: features unassigned at their own rank (e.g. the
  1.x `...f__lachnospiraceae.other`) lump unrelated organisms and are dropped.
* `pool_effects(method = "REML", test = "adhoc")`, FDR within rank, part and term
  (see F7). `method = "DL", test = "z", adjust_by = c("part", "term")` reproduces 1.x.

## 4. Meta-analysis (sex difference, four studies)

Pooling the stored per-study results of the four studies (Bangladesh, Haiti, USA CA/FL,
USA NC) with `pool_effects(method = "DL", test = "z")` reproduces the stored 1.x pooled
result (`metab.sex$random`) for the sex term: same 185 taxa; estimate, se, confidence
limits, p and adjusted p agree to about 1e-12. The paper's conclusions reproduce: no
phylum (0 of 7) with pooled p < 0.05; four genera with p < 0.05 (Coprococcus,
Eubacterium, Megamonas, Pseudoramibacter); only Coprococcus after FDR.

**F6. 1.x pooled taxa from a single study.** `meta.taxa()` counted the rows a taxon has,
not the studies with an estimate for the term, so for `bfNo_BF` 82 of 185 "pooled" taxa
came from one study and 4 had no estimate. The sex term is unaffected (all four studies
report it). 2.0 counts studies with an estimate.

**F7. Plain Hartung-Knapp is unsafe with very few studies.** With the 1.x data, plain
Hartung-Knapp (`test = "knha"`) made six genera significant instead of four: four of them
had only two studies with near-identical estimates, normal-theory p 0.22 to 0.86 but
Hartung-Knapp p 3e-5 to 7e-5. The ad hoc variant (`test = "adhoc"`, the 2.0 default) gave
0.43 to 0.89. With REML, the ad hoc test and FDR within the genus rank, only Coprococcus
remains (p 1.4e-4, FDR 0.012); the three weaker 1.x hits do not. The headline conclusion
holds; the secondary hits do not survive more appropriate small-sample inference.

**F8. Unassigned lineages as "taxa".** Several apparent genus-level hits were pseudo-taxa
such as `...o__clostridiales.other.other`, aggregates of unrelated organisms.

**F9. Only Bangladesh can be refitted from raw data.** The sample data (gender, feeding,
subject) of the Haiti, CA/FL and NC studies are not in the repository (see
`dev/LEGACY-DATA.md`). A hybrid end-to-end check refits Bangladesh and pools it with the stored
results of the other three.

**F10. Unexplained standard-error differences.** Refitting the Bangladesh sex model
(`~ gender + bf + age.sample`, longitudinal, 20 iterations, as documented) reproduces the
stored estimates to 1e-13 for all 124 taxa, but not the standard errors: they differ in 77
taxa (median ratio 1.02, up to about 18% larger in the refit; significance at 0.05 differs in
4 taxa). This is not due to the GAMLSS version (5.0-5 and 5.5-5 give identical values) and not
to the standard-error method (`vcov` falls back to `qr` in these fits and neither matches).
The cause is unknown. The stored Bangladesh results of the other examples do reproduce exactly.
It does not change the meta-analysis conclusions: with the stored or the refitted Bangladesh
results, under the 1.x and the 2.0 settings, the same genera are found (only Coprococcus
survives FDR; 4 genera at p < 0.05 with the 1.x settings, 1 with the 2.0 defaults).
Scripts: `data-raw/refit-sex-bangladesh.R`, `data-raw/diagnose-se-discrepancy/`.

## 5. LMAS (arcsine square root linear model)

`fit_study(engine = "lm_asin")` equals `taxa.compare(propmed.rel = "lm", transform =
"asin.sqrt")` for the mixed and the linear model (tested), and matches the stored results.
Tables 3 and 4: all estimates, limits, p-values and FDR reproduce once the p-value error is
accounted for.

**F11. The published mixed-model p-values carry the `1.96 * pnorm` error.** Every stored
LMAS p-value for a mixed model equals 0.98 times the correct two-sided p-value (ratio 0.98
for all 124 taxa of Example 1). Published Table 3 values: 0.0088, 0.0292, 0.0668, 0.0103;
correct: 0.0089, 0.0298, 0.0682, 0.0105. No conclusion changes; this is the bug fixed in 1.3.

**F12. Errors in the printed Table 5 (LMAS columns).** Example 3a, Bacteroidetes: the table
prints +0.06 with interval (-0.12, 0.01); the computed and stored value is -0.056, so the sign
is wrong. Example 3b (EBF > 2 months): the Actinobacteria and Firmicutes rows are swapped
(the printed "Actinobacteria 0.00, p 0.9626" is the computed Firmicutes row), and the printed
FDR values (0.9989, 0.9707, 0.9707, 0.9191) cannot be obtained as Benjamini-Hochberg over
the four displayed or all eight phylum-level features (all give 0.9626). No conclusion
changes: nothing is significant in that block.

## 6. How to rerun

```
Rscript data-raw/validate-2019-tables.R           # phylum level, about 25 s
FULL=1 Rscript data-raw/validate-2019-tables.R    # plus all taxa of Example 1, about 2.5 min
Rscript data-raw/convergence-sensitivity.R        # main phyla, 20 vs 500 iterations
Rscript data-raw/convergence-all-taxa.R           # sample of non-converged taxa, about 3 min
Rscript data-raw/validate-2019-meta.R             # pooling step, a few seconds
Rscript data-raw/validate-2019-lmas.R             # LMAS columns, about 1 min
Rscript data-raw/refit-sex-bangladesh.R           # Bangladesh refit plus pooling, about 3 min
```
