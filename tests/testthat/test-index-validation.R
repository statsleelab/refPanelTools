panel <- make_test_panel()

index_lines <- function() readLines(gzfile(panel$index))

# Writes a BGZF index built from `lines` and returns its path.
with_index <- function(lines, name) {
  path <- file.path(panel$dir, name)
  write_bgzf(path, paste0(lines, "\n"), lines_per_block = 3L)
  path
}

# Every function that reads the index goes through the same parser, so run the
# malformed cases against all of them.
index_readers <- list(
  extract_chr_data = function(idx, out)
    extract_chr_data(1, panel$num_pops, idx, panel$geno, out),
  extract_reg_data = function(idx, out)
    extract_reg_data(1, 0, 1e9, panel$num_pops, idx, panel$geno, out),
  extract_chr_pop_data = function(idx, out)
    extract_chr_pop_data(1, c("POPA"), idx, panel$geno, panel$pop_desc, out),
  extract_all_af1 = function(idx, out)
    extract_all_af1(1, idx, panel$geno, panel$pop_desc, out),
  simulate_af1_z = function(idx, out)
    simulate_af1_z(1, c("POPA"), 10L, idx, panel$geno, panel$pop_desc, out)
)

test_that("an fpos in scientific notation is rejected", {
  # fwrite() writes a large offset this way when the column is a double rather
  # than integer64. "1.2e+15" used to parse as 1, seeking a byte into the file
  # and dropping the first subject from every genotype string.
  lines <- index_lines()
  lines[1] <- sub(" [0-9]+$", " 1.2e+15", lines[1])
  idx <- with_index(lines, "idx_scientific.gz")

  for (nm in names(index_readers))
    expect_error(index_readers[[nm]](idx, tempfile()),
                 "unexpected field", info = nm)
})

test_that("a truncated index line is rejected", {
  # A short line left fpos at 0, so the SNP at offset 0 was silently emitted in
  # place of the intended one.
  lines <- index_lines()
  lines[2] <- "rs2 1 1100"
  idx <- with_index(lines, "idx_short.gz")

  for (nm in names(index_readers))
    expect_error(index_readers[[nm]](idx, tempfile()),
                 "malformed index line", info = nm)
})

test_that("a negative file offset is rejected", {
  lines <- index_lines()
  lines[1] <- sub(" [0-9]+$", " -5", lines[1])
  idx <- with_index(lines, "idx_negative.gz")
  expect_error(index_readers$extract_chr_data(idx, tempfile()),
               "negative file offset")
})

test_that("non-numeric chromosomes are skipped, not rejected", {
  # Panels may carry X, Y or MT rows. These functions select by chromosome
  # number and cannot address such rows, but their presence must not stop a run.
  idx <- with_index(c(index_lines(), "rsX X 9000 A G 0.10000 0"), "idx_chrX.gz")
  out <- tempfile()
  expect_silent(extract_chr_data(1, panel$num_pops, idx, panel$geno, out))
  expect_identical(readLines(out), panel$geno_lines[panel$chr == 1L])
})

test_that("a valid index still reads cleanly", {
  idx <- with_index(index_lines(), "idx_roundtrip.gz")
  out <- tempfile()
  extract_chr_data(2, panel$num_pops, idx, panel$geno, out)
  expect_identical(readLines(out), panel$geno_lines[panel$chr == 2L])
})

test_that("an unwritable output path is an error", {
  bad_out <- file.path(panel$dir, "no_such_directory", "out.txt")
  expect_error(cal_af1ref(panel$geno, panel$num_pops, bad_out), "can't open output file")
  expect_error(indexer(panel$geno, bad_out), "can't open output file")
  for (nm in names(index_readers))
    expect_error(index_readers[[nm]](panel$index, bad_out),
                 "can't open output file", info = nm)
})

test_that("a plain gzip file reports that bgzip compression is required", {
  plain <- file.path(panel$dir, "plain_gzip.gz")
  con <- gzfile(plain, "w"); writeLines(index_lines(), con); close(con)
  expect_error(extract_chr_data(1, panel$num_pops, plain, panel$geno, tempfile()),
               "compressed with bgzip")
})

test_that("a failed read does not leak the BGZF handles", {
  # bgzf_close() was only reached by falling off the end of a function, so
  # every Rcpp::stop() in between leaked both open handles. Before the reader
  # became scope-owned, 200 failed calls leaked 400 descriptors and a session
  # doing this in a loop would eventually run out.
  skip_on_os("windows")
  skip_if_not(dir.exists("/dev/fd"))

  idx <- with_index(replace(index_lines(), 2, "rs2 1 1100"), "idx_leak.gz")
  open_fds <- function() length(list.files("/dev/fd", all.files = TRUE))

  before <- open_fds()
  for (i in seq_len(100))
    try(extract_chr_data(1, panel$num_pops, idx, panel$geno, tempfile()),
        silent = TRUE)
  expect_lt(open_fds() - before, 10L)
})
