# Validation of the 2.0 engine against Ho et al. (2019)

Status: first gate of the 2.0 plan. The 2.0 `fit_study()` engine `bezi_gamlss`
(the 1.x GAMLSS-BEZI model) is compared with (a) the result objects stored in
the original GitHub repository (`data-raw/legacy-1.1/results/`) and (b) the numbers
printed in Tables 3 to 5 of the paper. Scripts: `data-raw/validate-2019-tables.R`,
`data-raw/convergence-sensitivity.R`, `data-raw/convergence-all-taxa.R`.
Environment: R 4.3.3, gamlss 5.5-5 (the paper used gamlss 5.0-5).

## 1. Reproduction

| Comparison | Result |
|---|---|
| New engine vs numbers printed in Tables 3, 4, 5 (GAMLSS-BEZI columns): estimate, both 95% limits, p, FDR for 16 phylum rows | **80 of 80 match** (to the printed 2 or 4 decimals) |
| New engine vs stored 1.x results, phylum level, 4 examples | max absolute difference about 1e-15 |
| New engine vs stored 1.x results, all 124 taxa of Example 1 | same 124 taxa; max absolute difference 2e-11 |
| New engine vs 1.x `taxa.compare()`, same data (unit test) | equal to 1e-8 |

Not yet reproduced: the LMAS (arcsin-sqrt linear model) columns, the simulation
studies, and the four-study meta-analysis (Fig. 5).

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

## 3. How to rerun

```
Rscript data-raw/validate-2019-tables.R           # phylum level, about 25 s
FULL=1 Rscript data-raw/validate-2019-tables.R    # plus all taxa of Example 1, about 2.5 min
Rscript data-raw/convergence-sensitivity.R        # main phyla, 20 vs 500 iterations
Rscript data-raw/convergence-all-taxa.R           # sample of non-converged taxa, about 3 min
```
