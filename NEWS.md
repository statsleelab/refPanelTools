# refPanelTools 1.1.0

## New functions

* `read_region()`: A convenience wrapper around `extract_reg_data()` that
  returns extracted SNP data directly as a data.frame, combining genotype
  strings and per-population allele frequencies with SNP metadata (rsid, chr,
  bp, a1, a2) from the index file. Eliminates the manual parsing step
  previously required after calling `extract_reg_data()`.

## Improvements

* Added roxygen2 documentation for all exported functions. Help pages are now
  available via `?function_name` (e.g. `?extract_reg_data`,
  `?simulate_af1_z`).

* Eliminated ~250 lines of duplicated code between `simulate_af1_z()` and
  `simulate_af1_z_allchr()` by extracting a shared internal helper. External
  behaviour and function signatures are unchanged.

* Unified error handling across all C++ functions: replaced `std::cout` /
  `exit()` calls with `Rcpp::stop()` and `Rcpp::Rcout()` so that errors are
  propagated cleanly to R and can be caught with `tryCatch()`.

* Removed `stderr` usage from `bgzf.c` for CRAN compatibility.

## Package metadata

* Updated `DESCRIPTION`: added `URL`, `BugReports`, proper `Authors@R` field,
  `Suggests: testthat`, and bumped version to 1.1.0.

* `NAMESPACE` is now managed by roxygen2.

---

# refPanelTools 1.0

* Initial release.
