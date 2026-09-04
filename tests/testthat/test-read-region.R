panel <- make_test_panel()

test_that("read_region() joins SNP metadata onto the genotype data", {
  got <- read_region(1, 900, 1250, panel$index, panel$geno,
                     pop_desc_file = panel$pop_desc)

  expect_identical(got$rsid, panel$rsid[panel$chr == 1L])
  expect_identical(got$bp,   panel$bp[panel$chr == 1L])
  expect_equal(got$af1ref,   round(panel$af1ref[panel$chr == 1L], 5))
  expect_identical(got$a1,   rep("A", 3))
  expect_identical(got$a2,   rep("G", 3))
})

test_that("read_region() names genotype and AF1 columns after the populations", {
  got <- read_region(1, 900, 1250, panel$index, panel$geno,
                     pop_desc_file = panel$pop_desc)
  expect_named(got, c("rsid", "chr", "bp", "a1", "a2", "af1ref",
                      "geno_POPA", "geno_POPB", "geno_POPC",
                      "af1_POPA", "af1_POPB", "af1_POPC"))
  expect_identical(got$geno_POPA, c("0120", "1111", "0000"))
  expect_identical(got$geno_POPC, c("012012", "111111", "000000"))
})

test_that("read_region() preserves leading zeros in genotype strings", {
  # Genotypes are packed digit strings, so they must not be read as numbers.
  got <- read_region(1, 900, 1250, panel$index, panel$geno,
                     pop_desc_file = panel$pop_desc)
  expect_type(got$geno_POPA, "character")
  expect_identical(nchar(got$geno_POPA), rep(4L, 3))
  expect_identical(nchar(got$geno_POPC), rep(6L, 3))
})

test_that("read_region() falls back to generic names when given num_pops", {
  got <- read_region(2, 1250, 1450, panel$index, panel$geno,
                     num_pops = panel$num_pops)
  expect_true(all(c("geno_pop1", "af1_pop3") %in% names(got)))
  expect_identical(got$rsid, c("rs4", "rs5"))
})

test_that("read_region() needs either a description file or num_pops", {
  expect_error(read_region(1, 900, 1250, panel$index, panel$geno),
               "pop_desc_file.*num_pops")
})

test_that("read_region() returns no rows for an empty region", {
  got <- read_region(1, 1, 2, panel$index, panel$geno,
                     pop_desc_file = panel$pop_desc)
  expect_identical(nrow(got), 0L)
})

test_that("read_region() agrees with extract_reg_data()", {
  out <- tempfile()
  extract_reg_data(1, 900, 1250, panel$num_pops, panel$index, panel$geno, out)
  raw <- readLines(out)
  got <- read_region(1, 900, 1250, panel$index, panel$geno,
                     pop_desc_file = panel$pop_desc)
  expect_identical(nrow(got), length(raw))
  expect_identical(got$geno_POPA, vapply(strsplit(raw, " "), `[`, character(1), 1L))
})
