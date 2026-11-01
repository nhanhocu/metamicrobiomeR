## Status

metamicrobiomeR 2.0 is in development and has not been submitted to CRAN.

## Checks so far (not a submission record)

* local: Ubuntu 24.04, R 4.3.3, `R CMD check --as-cran` (without the CRAN incoming checks): no errors; the only
  warning is a locale message from the check environment, and one note about file timestamps from the sandbox clock
* tests and vignettes: all pass
* <!-- TODO before submitting: win-builder (devel and release), R-hub, GitHub Actions (untested workflow) -->
* <!-- TODO: decide CRAN or Bioconductor (the importers use SummarizedExperiment and phyloseq) -->

## Downstream dependencies

None known.
