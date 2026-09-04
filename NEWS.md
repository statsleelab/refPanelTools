# refPanelTools 1.1.0

## Bug fixes

* `simulate_af1_z()` / `simulate_af1_z_allchr()`: **`num_sim_vec` was matched to
  populations in the order they appear in the population description file, not
  in the order given in `pop_vec`.** Whenever `pop_vec` was supplied in a
  different order, the per-population sample sizes were silently permuted and
  the simulated AF1 and Z-scores were wrong. `num_sim_vec[i]` is now always the
  sample size for `pop_vec[i]`. Results are unchanged for callers who already
  listed `pop_vec` in population-description-file order.

* `simulate_af1_z()` / `simulate_af1_z_allchr()`: added input validation --
  duplicate entries in `pop_vec`, non-positive values in `num_sim_vec`, an
  empty `pop_vec`, and a total sample size below 3 (the minimum needed for the
  Z-score denominator) are now errors instead of silently producing `nan`.

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
