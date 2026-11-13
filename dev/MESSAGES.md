# Messages: style guide and audit

## Principles

The package must be quiet when nothing notable happens, and specific when something does.

| Kind | Use it for | Never for |
|---|---|---|
| **Error** (`stop`) | The request cannot be carried out | Anything the function can sensibly do instead |
| **Warning** (`warning`) | Something likely to affect the result: failed or non-converged fits, estimates not used, unsafe pooling, no feature left to analyse | Routine information |
| **Message** (`.inform`) | A decision made for the user, or data silently left out: samples with missing values, features below the filters, features not pooled, unassigned reads, a subject effect added by default | Anything the user asked for explicitly |
| **Printed output** (`print` methods) | Show what an object is: a header, then the first rows | Raw data-frame dumps without context |

## Wording rules

1. Say **what is wrong, what was received, and what to do**, in at most two sentences.
2. **Lead with the count**, then the thing: "6 of 40 samples have missing values", not "Some samples...".
3. **Name the specific IDs**, at most three, then "(n more)": `.show()`.
4. Argument names in single quotes (`'samples'`), functions with parentheses (`agglomerate()`), values in the data in plain text.
5. Proper singular/plural (`.plural()`, `.word()`); no "(s)". "an mm_study" (read "em-em").
6. Sentence case, ending with a full stop. No jargon without a pointer (`?pool_effects`, `see attr(, "fit_info")`).
7. Warnings and errors use `call. = FALSE`; the message says which function and argument.

## Controlling messages

* `options(metamicrobiomeR.quiet = TRUE)` silences messages (not warnings); `suppressMessages()` works too.
* A message that would repeat on every call (the default subject effect, normalised function tables,
  the two-study Hartung-Knapp note) is shown **once per session**.
* Third-party console output (lme4 singular-fit notes, glmmTMB and GAMLSS messages) is captured per
  feature and reported once in aggregate; the details are in `attr(effects, "fit_info")` (`first_warning`, `singular`).

## Helpers (`R/01-messages.R`)

`.inform()`, `.show()`, `.plural()`, `.word()`, `.cap()`, `.pct()`, `.check_class()`, `.check_list_of()`,
`.id_mismatch()`, `.missing_vars()`, `.need()` (names the install command, `BiocManager::install()` for Bioconductor
packages), `.check_comparable()` (one engine and one scale per part and term, used by `bind_effects()`,
`pool_effects()` and `meta_regress()`).

## Catalogue of informational messages and warnings

| Where | When | Text (abridged) |
|---|---|---|
| `fit_study()` | features filtered | Fitting 5 of 8 features (2 unassigned, 1 below the prevalence and abundance filters). |
| | no feature passes (warning) | No feature passed the filters (present in at least 5% ...); lower 'min_prevalence' or 'min_mean'. |
| | missing covariates | 6 of 40 samples have missing values in the model variables and are left out. |
| | subject default (once) | Using 'person' as a subject random intercept (set longitudinal = FALSE to ignore it). |
| | singular mixed fits | 5 of 8 mixed-model fits were singular (the subject variance is estimated as zero), so ... |
| | failures, no convergence (warning) | 3 of 8 fits failed (first error: ...); see attr(, "fit_info"). |
| `pool_effects()` | features not pooled | 2 features of 12 not pooled (estimated in fewer than 2 studies or ...). |
| | two studies (once) | 1 feature was pooled from only 2 studies: the Hartung-Knapp test then has 1 degree of freedom ... |
| | estimates unusable (warning) | 2 estimates not used: 1 with a missing estimate or standard error, 1 with ... |
| `new_effect_table()` | zero standard errors (warning) | 1 standard error is exactly zero; such estimates cannot be weighted ... |
| `as_study()` | zero-depth samples (warning) | 1 sample has zero total abundance (S2); kept, but left out of models. |
| | taxonomy without ranks (warning) | 'taxonomy' has no rank columns (kingdom, ..., species; any case, 'domain' counts as kingdom) ... |
| `agglomerate()` | unassigned features | 2 of 5 features (14% of the abundance) have no assignment at 'genus'; kept as 'unassigned_genus' ... |
| `merge_studies()` | summary; low overlap (warning) | Merged 3 studies: 120 features (of 200, those present in all studies) in 330 samples. |
| `alpha_diversity()` | rarefaction | Rarefying to 25,913 reads: 30 of 40 samples with fewer reads are dropped. |
| `fit_beta()` | samples left out | 3 of 20 samples are left out (missing model variables or no abundance). |
| `adjust_hierarchical()` | testing summary | 5 of 12 features were tested; 7 were not because their parent was not rejected. ... |
| `meta_regress()` | skipped features | 1 of 2 features skipped (fewer studies than model coefficients, or ...). |
| `filter_effects()` | what was dropped | Dropped 1 of 6 estimates before pooling: 1 with an implausibly small standard error ... |
| `read_dada2()` | order assumed (warning) | The row names of 'taxa' do not match the sequences of 'seqtab'; assuming the rows are in the same order. |
| `read_qiime2()` | no taxonomy | 2 of 3 features have no row in the taxonomy table and are left unassigned. |
| `read_metaphlan()` | level kept | Kept 1 of 2 rows: the clades at level 'species'. |
| `read_function_table()` | normalised (once) | Each sample is normalised to proportions, so functions are analysed as relative abundances. |

## Audit (14 scenarios, run before and after)

| Scenario | Before | After |
|---|---|---|
| Mixed model with a near-zero subject variance | 5 identical lme4 messages, one per feature | the default subject effect noted once, singular fits reported once |
| No feature passes the filters | silent empty table | warning that says why and what to lower |
| 6 of 40 samples with missing values | silent | message with the counts |
| Subject column set | random effect added silently | said once per session |
| Pooling a table with two engines | ran | error naming the part, term and engines |
| Taxonomy column called `Domain` | silently ignored | recognised as the kingdom |
| QIIME 2 features without taxonomy; rarefaction dropping samples | silent | message with the counts |
| Sample IDs that differ | "Row names of 'samples' must match..." | names the IDs that differ and in which table |
| Wrong object passed | "'study' must be an mm_study object." | adds where to get one and what was received |
| Printing an effect table | 18 raw lines | a header, then six rows |

Every message above has a test in `tests/testthat/test-01-messages.R`, including that nothing is printed when
nothing notable happens. Add a test with any new message.
