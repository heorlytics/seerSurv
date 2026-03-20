# CRAN submission comments — seerSurv 0.1.0

## Test environments

* local: macOS 15 (Sequoia), R 4.4.2
* GitHub Actions: ubuntu-latest (R release, devel, oldrel-1), windows-latest (release), macos-latest (release)
* win-builder: R-devel
* rhub: Fedora Linux (R-devel), Windows Server 2022 (R-devel)

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs on any test environment.

## Downstream dependencies

This is a new submission; there are no downstream packages.

## Comments to CRAN reviewers

* The package bundles a CSV-derived `.rda` dataset (`lifetable_seer`) and a
  hard-coded tibble (`tumour_data_seer`).  Both are sourced from the US-CDC
  National Vital Statistics system and the NCI SEER-17 database, which are
  both in the public domain.

* Two long-running integration tests are wrapped in `skip_on_cran()` to keep
  the `R CMD check --as-cran` time well under 5 minutes.

* The examples in `run_tumour_analysis()` and `blend_survival()` are wrapped
  in `\donttest{}` because they require `flexsurv` to fit multiple models and
  take ~5 seconds each.
