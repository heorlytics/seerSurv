# Contributing to seerSurv

We welcome contributions of all kinds — bug reports, documentation
improvements, new features, and additional cancer-type data.

---

## Reporting bugs

Please open an issue on
[GitHub](https://github.com/pharmacoevidence/seerSurv/issues) and include:

1. A minimal reproducible example (use `reprex::reprex()`).
2. The output of `sessionInfo()`.
3. The exact error message or unexpected output.

---

## Proposing changes

1. Fork the repository and create a branch from `main`.
2. Install development dependencies:

   ```r
   install.packages(c("devtools", "roxygen2", "testthat", "covr", "lintr"))
   ```

3. Make your changes.  Follow the [tidyverse style guide](https://style.tidyverse.org/).
4. Document new functions with **roxygen2** tags and run `devtools::document()`.
5. Add or update **testthat** tests.  Aim for > 90 % line coverage.
6. Run the full check suite and confirm it passes with no warnings or notes:

   ```r
   devtools::check()
   ```

7. Submit a pull request with a clear description of what changed and why.

---

## Code style

* Indentation: 2 spaces (no tabs).
* Line length: ≤ 100 characters.
* Avoid `T` / `F` — use `TRUE` / `FALSE`.
* Use `\(x)` lambda syntax (R ≥ 4.1) rather than `function(x)` for
  short anonymous functions.
* All user-facing functions must have complete roxygen2 documentation
  including at least one `@examples` block.
* Internal helpers (prefixed with `.`) do not require exported documentation
  but should have inline comments.

---

## Testing

```r
devtools::test()
covr::package_coverage()
```

---

## Code of conduct

This project follows the
[Contributor Covenant Code of Conduct v2.1](https://www.contributor-covenant.org/version/2/1/code_of_conduct/).
By participating, you agree to abide by its terms.
