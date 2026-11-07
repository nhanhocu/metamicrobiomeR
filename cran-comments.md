## Submission

This is a resubmission of metamicrobiomeR, which was archived on 2026-03-25 because the issues reported for version 1.2
were not corrected in time. Version 2.0.0 corrects them: the global variable `id` in `metatab.show()` is declared, and
`dplyr` is imported by name rather than as a whole. It also adds new functions, while the functions of version 1.x stay
available with the same arguments.

The maintainer changes from Nhan Ho (nhanhocumc@gmail.com) to Payton Yau (tungon@gmail.com). Nhan Ho remains an author and
will confirm the change by email.

## Test environments

* GitHub Actions: ubuntu-latest (R-devel, release, oldrel-1), windows-latest (release), macos-latest (release)
* <!-- TODO before submitting: win-builder (R-devel and release), and R-hub -->

## R CMD check results

0 errors | 0 warnings | 0 notes on GitHub Actions, apart from the notes that CRAN adds on submission:

* New maintainer, and the package was archived on CRAN.
* <!-- TODO: add any other note from win-builder, with a one-line explanation -->

## Downstream dependencies

There are no reverse dependencies on CRAN.
