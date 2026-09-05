# refPanelTools 1.1.0

## Testing

* Added a `testthat` suite (75 assertions) covering every exported function.
  The BGZF fixtures are generated in `tests/testthat/helper-panel.R` rather
  than checked in as binaries, so the panel format the C++ code expects is
  written out in readable form. Each of the bugs fixed in this release has a
  regression test, and those tests were checked against the unfixed code to
  confirm they actually fail on it.

## Documentation

* `indexer()`'s documentation described it as producing the index the other
  functions read. It does not: it writes two columns (`fpos line_length`), and
  the extract and simulate functions read seven
  (`rsid chr bp a1 a2 af1ref fpos`). Its actual role is to recompute offsets
  after a genotype file has been derived from an existing panel, so that the
  `fpos` column can be replaced while the SNP metadata is carried over. The
  help page now says so and works the per-chromosome split through as an
  example. `extract_chr_data()` no longer claims its index is `indexer()`'s
  output, and every function that reads an index now names the seven columns.

## Bug fixes

* `read_region()` read the entire index file with `read.table()` on every
  call, so the cost tracked the size of the index rather than the size of the
  region: a 40-SNP window cost the same as a 40,000-SNP one. The metadata it
  was recovering had already been parsed on the C++ side and thrown away. A new
  internal `read_reg_records()` keeps those fields as it parses them and
  returns them with the genotypes in one pass, and `read_region()` no longer
  touches the index from R. Measured 3.4x faster on an 800,000-row index, with
  the memory for the index copy gone.

* `read_region()` returned the `af1_*` columns as character, because the
  `colClasses = "character"` needed to protect genotype strings from losing
  their leading zeros was applied to every column. They are numeric now;
  genotype strings are still character.

* `read_region()` returned a bare `data.frame()` with no columns when a region
  held no SNPs, so `$rsid` and `rbind()` across regions broke on the empty
  case. It now returns zero rows with the same columns as a non-empty result.

* `read_region()` now rejects a `num_pops` that does not match the panel --
  which would otherwise shift every genotype and AF1 column -- and a reversed
  base-pair range.

* The index file was parsed with unchecked `>>` extractions, turning two
  realistic ways of corrupting an index into wrong answers rather than errors.
  A truncated line left `fpos` at 0, so the SNP stored at offset 0 was emitted
  in place of the intended one. An `fpos` in scientific notation -- what
  `fwrite()` produces for a large offset when \pkg{bit64} is not installed to
  keep the column `integer64` -- parsed as its leading digit, seeking a few
  bytes into the file and shifting every genotype string by one subject.
  Parsing now goes through one checked helper used by all five readers.
  Chromosomes that are not numbers (X, Y, MT) are still skipped rather than
  rejected, as before.

* Output files were opened without checking that the open succeeded, so an
  unwritable path looked like a successful run that produced nothing. All
  seven call sites now report the failure.

* Reading a plain gzip file failed with "can't read 1-th character from BGZF
  file". The message now says the file must be compressed with `bgzip`.

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

* `simulate_af1_z()` / `simulate_af1_z_allchr()`: the bootstrap sample and the
  null phenotype were drawn from `std::mt19937` seeded by `std::random_device`,
  so `set.seed()` had no effect and no run could be reproduced -- a simulation
  whose output could not be regenerated or checked. Both now draw from R's own
  generator (`norm_rand()` and `R_unif_index()`), so a run is reproducible with
  `set.seed()`, including across sessions. Runs without an explicit seed still
  differ from each other, as before. `R_unif_index()` is the draw `sample()`
  uses and avoids the modulo bias of scaling a uniform deviate by hand.

* `simulate_af1_z()` / `simulate_af1_z_allchr()`: a population declared with
  zero subjects in the population description file is now an error rather than
  an empty sampling range.

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
