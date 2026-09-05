panel <- make_test_panel()

test_that("test_gz_file() prints the first line of a BGZF file", {
  expect_output(test_gz_file(panel$geno), fixed = TRUE,
                regexp = panel$geno_lines[1])
})

test_that("opening a missing file is an error", {
  missing <- file.path(panel$dir, "does_not_exist.gz")
  expect_error(test_gz_file(missing), "can't open")
  expect_error(cal_af1ref(missing, panel$num_pops, tempfile()), "can't open")
  expect_error(extract_chr_data(1, panel$num_pops, missing, panel$geno, tempfile()),
               "can't open")
})

test_that("cal_af1ref() computes AF1 across all populations", {
  out <- tempfile()
  cal_af1ref(panel$geno, panel$num_pops, out)
  expect_equal(as.numeric(readLines(out)), round(panel$af1ref, 5))
})

test_that("cal_af1ref() rounds to nearest, not up", {
  # rs1 has 13 alternate alleles over 15 subjects: 0.4333333..., which rounds
  # to 0.43333. std::ceil() used to report 0.43334.
  out <- tempfile()
  cal_af1ref(panel$geno, panel$num_pops, out)
  expect_identical(readLines(out)[1], "0.43333")
})

test_that("extract_chr_data() returns every SNP on the chromosome, in order", {
  for (chr in c(1L, 2L)) {
    out <- tempfile()
    extract_chr_data(chr, panel$num_pops, panel$index, panel$geno, out)
    expect_identical(readLines(out), panel$geno_lines[panel$chr == chr])
  }
})

test_that("extract_chr_data() seeks correctly across BGZF block boundaries", {
  # The fixture packs two SNPs per block, so chromosome 2 starts mid-file and
  # its two SNPs sit in different blocks. A wrong virtual offset would surface
  # here rather than on chromosome 1.
  out <- tempfile()
  extract_chr_data(2, panel$num_pops, panel$index, panel$geno, out)
  expect_identical(readLines(out), panel$geno_lines[panel$chr == 2L])
  expect_gt(sum(panel$chr == 2L), 1L)   # the fixture must span blocks
})

test_that("extract_chr_data() writes nothing for an absent chromosome", {
  out <- tempfile()
  extract_chr_data(99, panel$num_pops, panel$index, panel$geno, out)
  expect_identical(readLines(out), character(0))
})

test_that("extract_reg_data() filters on an inclusive bp range", {
  out <- tempfile()
  extract_reg_data(1, 1000, 1100, panel$num_pops, panel$index, panel$geno, out)
  expect_identical(readLines(out), panel$geno_lines[1:2])

  extract_reg_data(1, 1001, 1099, panel$num_pops, panel$index, panel$geno, out)
  expect_identical(readLines(out), character(0))
})

test_that("extract_chr_pop_data() keeps only the requested populations", {
  out <- tempfile()
  extract_chr_pop_data(1, c("POPA", "POPC"), panel$index, panel$geno,
                       panel$pop_desc, out)
  got <- strsplit(trimws(readLines(out)), " +")
  # two genotype strings then two AF1 values, in population-description order
  expect_identical(vapply(got, `[`, character(1), 1L),
                   c("0120", "1111", "0000"))
  expect_identical(vapply(got, `[`, character(1), 2L),
                   c("012012", "111111", "000000"))
  expect_identical(vapply(got, length, integer(1)), rep(4L, 3))
})

test_that("extract_chr_pop_data() matches population names case-insensitively", {
  upper <- tempfile(); lower <- tempfile()
  extract_chr_pop_data(1, c("POPA"), panel$index, panel$geno, panel$pop_desc, upper)
  extract_chr_pop_data(1, c("popa"), panel$index, panel$geno, panel$pop_desc, lower)
  expect_identical(readLines(upper), readLines(lower))
})

test_that("extract_all_af1() writes per-population AF1 in panel order", {
  out <- tempfile()
  extract_all_af1(1, panel$index, panel$geno, panel$pop_desc, out)
  got <- do.call(rbind, lapply(strsplit(trimws(readLines(out)), " +"), as.numeric))
  expect_equal(got[1, ], c(0.375, 0.4, 0.5))
  expect_equal(got[3, ], c(0, 0, 0))
  expect_identical(dim(got), c(3L, 3L))
})

test_that("indexer() records a byte offset and length for every line", {
  out <- tempfile()
  indexer(panel$geno, out)
  got <- read.table(out, col.names = c("fpos", "len"))
  expect_identical(nrow(got), length(panel$geno_lines))
  expect_identical(got$len, nchar(panel$geno_lines))
  # Note: this two-column layout is *not* the seven-column index the extract_*
  # and simulate_* functions consume (rsid chr bp a1 a2 af1ref fpos).
  expect_identical(ncol(got), 2L)
})

test_that("a num_pops that does not match the panel is an error", {
  # num_pops used to be accepted and ignored by the extraction functions, so a
  # wrong value went unnoticed until the output was parsed downstream.
  wrong <- panel$num_pops + 1L
  expect_error(extract_chr_data(1, wrong, panel$index, panel$geno, tempfile()),
               "genotype fields, expected")
  expect_error(extract_reg_data(1, 900, 1250, wrong, panel$index, panel$geno,
                                tempfile()),
               "genotype fields, expected")
  expect_error(cal_af1ref(panel$geno, wrong, tempfile()),
               "genotype fields, expected")

  # the right value still works
  expect_silent(extract_chr_data(1, panel$num_pops, panel$index, panel$geno,
                                 tempfile()))
})
