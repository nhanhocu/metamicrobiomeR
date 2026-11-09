## Submission

This is a resubmission of metamicrobiomeR, which was archived on 2026-03-25 because the issues reported for version 1.2
were not corrected in time. Version 2.0.0 corrects them: the global variable `id` in `metatab.show()` is declared, and
`dplyr` is imported by name rather than as a whole. We also fixed related problems found while preparing this release:

* `class(x) == ...` comparisons are replaced by `inherits()`.
* The unused `lmerTest` import is removed, and `gdata` and `matrixStats` are replaced by base R.
* Examples no longer download files.
* `meta.taxa()` sets its between-study estimator explicitly, because the default of the `meta` package changed and stopped
  the function on its own example data.

The release also adds new functions, while the functions of version 1.x stay available with the same arguments. All
changes are listed in NEWS.md.

The maintainer changes from Nhan Ho (nhanhocumc@gmail.com) to Payton Yau (tungon@gmail.com). Nhan Ho remains an author and
will confirm the change by email.

## Test environments

* Local: Linux, R 4.5, `R CMD check --as-cran`
* GitHub Actions: ubuntu-latest (R-devel, release, oldrel-1), windows-latest (release), macos-latest (release)
* R-hub: <!-- TODO before submitting: run from the Actions tab (linux, windows, macos) and note the result -->
* win-builder: <!-- TODO before submitting: devtools::check_win_devel() and check_win_release(); results go to the maintainer -->

## R CMD check results

0 errors | 0 warnings | 0 notes on GitHub Actions and locally, apart from the notes that CRAN adds on submission:

* New maintainer, and the package was archived on CRAN.
* <!-- TODO: add any other note from win-builder or R-hub, with a one-line explanation -->

## Downstream dependencies

There are no reverse dependencies on CRAN.
