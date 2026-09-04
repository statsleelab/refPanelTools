panel <- make_test_panel()

# Reads a simulation result back as a data.frame.
sim <- function(..., allchr = FALSE) {
  out <- tempfile()
  if (allchr) simulate_af1_z_allchr(..., ref_out_file = out)
  else        simulate_af1_z(..., ref_out_file = out)
  utils::read.table(out, header = TRUE, stringsAsFactors = FALSE)
}

# --- num_sim_vec is paired with pop_vec positionally -----------------------

test_that("num_sim_vec follows pop_vec order, not the description file's", {
  # At the discriminating SNP POPA is fixed at genotype 2 and POPC at 0, so the
  # simulated AF1 is exactly n_POPA / (n_POPA + n_POPC) -- it reports the
  # realised mix. POPA precedes POPC in the description file, so a regression
  # that consumes num_sim_vec in file order swaps the first pair of calls.
  set.seed(1)
  af1 <- function(pv, nv) {
    got <- sim(2L, pv, nv, panel$index, panel$geno, panel$pop_desc)
    got[got$rsid == PANEL_DISCRIMINATING_SNP, "sim_af1"]
  }

  # sim_af1 is written to five decimal places, so compare against the rounded
  # expectation exactly rather than with a tolerance.
  expect_identical(af1(c("POPC", "POPA"), c(50L, 1L)), round(1 / 51, 5))
  expect_identical(af1(c("POPA", "POPC"), c(1L, 50L)), round(1 / 51, 5))
  expect_identical(af1(c("POPA", "POPC"), c(50L, 1L)), round(50 / 51, 5))
  expect_identical(af1(c("POPC", "POPA"), c(1L, 50L)), round(50 / 51, 5))
})

test_that("pop_vec order does not change the result when sizes track it", {
  set.seed(99); a <- sim(1L, c("POPA", "POPB", "POPC"), c(10L, 20L, 30L),
                         panel$index, panel$geno, panel$pop_desc)
  set.seed(99); b <- sim(1L, c("POPC", "POPB", "POPA"), c(30L, 20L, 10L),
                         panel$index, panel$geno, panel$pop_desc)
  expect_equal(a, b)
})

# --- input validation ------------------------------------------------------

test_that("invalid population and sample-size arguments are rejected", {
  bad <- function(pv, nv) sim(1L, pv, nv, panel$index, panel$geno, panel$pop_desc)

  expect_error(bad(c("POPA", "POPA"), c(10L, 10L)), "more than once")
  expect_error(bad(c("NOPE"),         c(10L)),      "not found in reference")
  expect_error(bad(c("POPA", "POPC"), c(10L)),      "same length")
  expect_error(bad(character(0),      integer(0)),  "must not be empty")
  expect_error(bad(c("POPA", "POPC"), c(0L, 10L)),  "positive values")
  expect_error(bad(c("POPA", "POPC"), c(-5L, 10L)), "positive values")
  expect_error(bad(c("POPA"),         c(2L)),       "at least 3")
})

test_that("a population size that disagrees with the genotype data is an error", {
  # The description file claims 500 subjects for POPA; the genotype strings
  # hold 4. Indexing past the end of the string used to run to completion and
  # report nonsense, including negative allele frequencies.
  bad_desc <- write_bad_pop_desc(panel$dir)
  expect_error(
    sim(1L, c("POPA"), c(20L), panel$index, panel$geno, bad_desc),
    "declared to have 500 subjects"
  )
})

test_that("a population declared with no subjects is an error", {
  zero <- write_bad_pop_desc(panel$dir, num_subj = c(0L, 5L, 6L))
  expect_error(sim(1L, c("POPA"), c(10L), panel$index, panel$geno, zero),
               "no subjects")
})

# --- undefined Z-scores ----------------------------------------------------

test_that("SNPs with no genotype variance report NA rather than nan", {
  # rs2 is all heterozygous and rs3 is monomorphic: sxx == 0 for both, leaving
  # the slope and its standard error undefined.
  out <- tempfile()
  set.seed(3)
  simulate_af1_z(1L, c("POPA", "POPB", "POPC"), c(10L, 10L, 10L),
                 panel$index, panel$geno, panel$pop_desc, out)

  # Assert on the file text, not the parsed frame: read.table() turns "nan"
  # into NaN, and is.na(NaN) is TRUE, so a parsed check would pass against the
  # old output too.
  raw <- readLines(out)
  expect_false(any(grepl("nan", raw, ignore.case = TRUE)))
  fields <- strsplit(raw[-1], " ")
  z_text <- setNames(vapply(fields, `[`, character(1), 7L),
                     vapply(fields, `[`, character(1), 1L))
  expect_identical(unname(z_text[c("rs2", "rs3")]), c("NA", "NA"))

  got <- utils::read.table(out, header = TRUE, stringsAsFactors = FALSE)
  expect_true(is.numeric(got$sim_z))
  expect_false(any(is.nan(got$sim_z)))
  expect_true(all(is.na(got$sim_z[got$rsid %in% c("rs2", "rs3")])))
  expect_false(is.na(got$sim_z[got$rsid == "rs1"]))
  expect_false(any(is.na(got$sim_af1)))   # AF1 stays defined
})

# --- reproducibility -------------------------------------------------------

test_that("set.seed() makes a run reproducible", {
  args <- list(1L, c("POPA", "POPB", "POPC"), c(10L, 10L, 10L),
               panel$index, panel$geno, panel$pop_desc)
  set.seed(42); a <- do.call(sim, args)
  set.seed(42); b <- do.call(sim, args)
  expect_equal(a, b)

  set.seed(43); c3 <- do.call(sim, args)
  expect_false(isTRUE(all.equal(a$sim_z, c3$sim_z)))
})

test_that("the RNG stream advances, so consecutive calls differ", {
  args <- list(1L, c("POPA", "POPB", "POPC"), c(10L, 10L, 10L),
               panel$index, panel$geno, panel$pop_desc)
  set.seed(7)
  a <- do.call(sim, args)
  b <- do.call(sim, args)
  expect_false(isTRUE(all.equal(a$sim_z, b$sim_z)))

  # and the draws are visible to R afterwards
  set.seed(7); before <- runif(1)
  set.seed(7); invisible(do.call(sim, args)); after <- runif(1)
  expect_false(isTRUE(all.equal(before, after)))
})

# --- output shape and statistical behaviour --------------------------------

test_that("simulate_af1_z() writes the documented columns for one chromosome", {
  set.seed(5)
  got <- sim(2L, c("POPA"), c(10L), panel$index, panel$geno, panel$pop_desc)
  expect_named(got, c("rsid", "chr", "bp", "a1", "a2", "sim_af1", "sim_z"))
  expect_identical(got$rsid, panel$rsid[panel$chr == 2L])
  expect_identical(got$bp,   panel$bp[panel$chr == 2L])
  expect_true(all(got$chr == 2L))
})

test_that("simulate_af1_z_allchr() covers every chromosome", {
  set.seed(5)
  got <- sim(c("POPA"), c(10L), panel$index, panel$geno, panel$pop_desc,
             allchr = TRUE)
  expect_identical(got$rsid, panel$rsid)
  expect_identical(got$chr,  panel$chr)
})

test_that("simulated AF1 converges on the population's true frequency", {
  # rs1 POPA is "0120": three alternate alleles over four subjects, AF1 0.375.
  # A sampler drawing outside [0, n) or skipping an index would miss this.
  set.seed(2024)
  got <- sim(1L, c("POPA"), c(50000L), panel$index, panel$geno, panel$pop_desc)
  expect_equal(got[1, "sim_af1"], 0.375, tolerance = 0.01)
})

test_that("null Z-scores are approximately standard normal", {
  # The phenotype is independent of genotype, so across independent runs the
  # Z-score for a given SNP should be N(0, 1). This exercises the whole
  # pipeline: sampling, the OLS fit, and the rounding.
  skip_on_cran()
  set.seed(11)
  z <- vapply(seq_len(200), function(i) {
    sim(1L, c("POPA", "POPB", "POPC"), c(70L, 70L, 60L),
        panel$index, panel$geno, panel$pop_desc)[1, "sim_z"]
  }, numeric(1))
  expect_lt(abs(mean(z)), 4 / sqrt(length(z)))   # ~4 SE
  expect_gt(sd(z), 0.8)
  expect_lt(sd(z), 1.2)
})
