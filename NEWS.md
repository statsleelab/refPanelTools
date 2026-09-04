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

* `bgzf.c`: `packInt16()`, `unpackInt16()` and `packInt32()` were declared
  `inline` without `static`, which under C99 and later emits no external
  definition. The package linked only because these were always inlined at
  `-O2`; unoptimised builds (including `devtools::load_all()` / `document()`)
  failed with `symbol not found: _packInt16`. Now declared `static inline`.

* `simulate_af1_z()` / `simulate_af1_z_allchr()`: the bootstrap sampling
  indices are drawn from `[0, num_subj)` using the subject count declared in
  the population description file, but the genotype string was indexed without
  checking its actual length. A declared count larger than the real one read
  past the end of the string -- undefined behaviour that silently produced
  nonsense (for example a negative `sim_af1`). The declared and actual subject
  counts are now required to agree, and a genotype record with too few fields
  is an error.

* `simulate_af1_z()` / `simulate_af1_z_allchr()`: SNPs that are monomorphic in
  the bootstrap sample have `sxx == 0`, leaving the slope and its standard
  error undefined; the literal string `nan` was written to `sim_z`. These are
  now reported as `NA`. Monomorphic SNPs are common in reference panels.

* `cal_af1ref()`, `simulate_af1_z()`, `simulate_af1_z_allchr()`: values were
  passed through `std::ceil(x * 1e5) / 1e5`, which rounds up rather than to
  nearest and shifts every value by +5e-06 on average. For example an AF1 of
  0.4333333 was reported as 0.43334 instead of 0.43333. For Z-scores the
  distortion is asymmetric, inflating positive values and shrinking negative
  ones. Now uses `std::round()`. **Output files produced by earlier versions
  differ in the fifth decimal place and should be regenerated.**

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
