# metamicrobiomeR

Differential abundance analysis and **meta-analysis across microbiome studies**.
Originally described in Ho et al. (2019), *BMC Bioinformatics* 20:188,
<https://doi.org/10.1186/s12859-019-2744-2>.

> **Status: version 2.0 is under development and is not on CRAN.** The package is
> being redesigned; the API will change until the first 2.0 release.

## What is changing in 2.0

Version 2.0 is a redesign and a clean break from 1.x. The input format of 1.x
(one wide table with QIIME 1 style column names) is replaced by a validated
study object with taxonomy in rank columns, optional counts and study-level
metadata, and by an effect table that connects per-study models to
meta-analysis. See [`dev/SPEC.md`](dev/SPEC.md) for the draft specification.

Planned: count-based models that use sequencing depth, rank-aware analysis,
taxonomy harmonisation across databases, heterogeneity-aware pooling with
moderators, and import from current formats (phyloseq, TreeSummarizedExperiment,
BIOM, DADA2, MetaPhlAn).

## Branches

| Branch | Content |
|---|---|
| `v2-dev` | Development of version 2.0 |
| `legacy-1.x` | The 1.x interface with bug fixes (version 1.3); the last release with the old API |
| `sync-cran-1.2` | CRAN 1.2 sources brought into this repository |

## Install

```r
# development version of 2.0
remotes::install_github("paytonyau/metamicrobiomeR@v2-dev")

# last release with the 1.x interface
remotes::install_github("paytonyau/metamicrobiomeR@legacy-1.x")
```

## Validation material

`data-raw/legacy-1.1/` holds the raw Bangladesh QIIME outputs and the 1.x result
objects from the original repository, used to check that 2.0 reproduces the 2019
results.

## Authors and licence

Nhan Ho (author), Payton Yau (maintainer). GPL (>= 2).
