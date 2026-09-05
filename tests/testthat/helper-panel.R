# Builds a small BGZF-compressed reference panel on disk so the extraction and
# simulation functions can be exercised end to end.
#
# The fixtures are generated rather than checked in: a binary blob in the
# repository would be opaque, and writing the format out here doubles as a
# specification of what the C++ code expects to read.

# --- BGZF writing ----------------------------------------------------------

# One BGZF block: an 18-byte gzip header carrying the "BC" extra field, the
# raw deflate stream, then CRC32 and ISIZE. R has no raw-deflate or CRC32
# primitive, so both are lifted out of a gzip stream written by gzfile().
bgzf_block <- function(payload) {
  # A BGZF block holds at most 64 KB either side of the compression: the reader
  # inflates into a fixed 64 KB buffer, and BSIZE is a 16-bit field. Exceeding
  # the first makes inflate fail, the second wraps silently; both surface as
  # "this file is not BGZF", which is a confusing way to learn the block is
  # simply too big.
  if (length(payload) > 65535L)
    stop("BGZF payload too large (", length(payload),
         " bytes uncompressed); lower lines_per_block")

  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  con <- gzfile(tmp, "wb")
  if (length(payload)) writeBin(payload, con)
  close(con)

  gz      <- readBin(tmp, "raw", file.info(tmp)$size)
  deflate <- gz[11:(length(gz) - 8)]   # strip the 10-byte header
  trailer <- gz[(length(gz) - 7):length(gz)]  # CRC32 + ISIZE, already LE

  bsize <- 18L + length(deflate) + 8L - 1L    # total block size minus one
  if (bsize > 65535L)
    stop("BGZF block too large (", bsize + 1L, " bytes compressed); ",
         "lower lines_per_block")
  header <- as.raw(c(
    0x1f, 0x8b, 0x08, 0x04,                   # magic, deflate, FEXTRA
    0x00, 0x00, 0x00, 0x00,                   # MTIME
    0x00, 0xff,                               # XFL, OS
    0x06, 0x00,                               # XLEN = 6
    0x42, 0x43,                               # subfield id "BC"
    0x02, 0x00,                               # SLEN = 2
    bsize %% 256L, bsize %/% 256L             # BSIZE, little endian
  ))
  c(header, deflate, trailer)
}

# Writes `lines` as a BGZF file and returns the virtual file offset of each
# line: the compressed offset of its block in the high 48 bits, its offset
# within the uncompressed block in the low 16. These are the `fpos` values the
# index file carries and bgzf_seek() consumes.
write_bgzf <- function(path, lines, lines_per_block = 2L) {
  stopifnot(length(lines) > 0L)
  out     <- raw(0)
  offsets <- numeric(length(lines))

  for (start in seq(1L, length(lines), by = lines_per_block)) {
    idx     <- start:min(start + lines_per_block - 1L, length(lines))
    coffset <- length(out)
    uoffset <- 0L
    for (i in idx) {
      offsets[i] <- coffset * 65536 + uoffset
      uoffset    <- uoffset + nchar(lines[i], type = "bytes")
    }
    out <- c(out, bgzf_block(charToRaw(paste0(lines[idx], collapse = ""))))
  }
  out <- c(out, bgzf_block(raw(0)))   # BGZF end-of-file marker block

  writeBin(out, path)
  offsets
}

# --- the test panel --------------------------------------------------------

# Three populations of different sizes, five SNPs over two chromosomes.
# Genotypes are chosen so that every expected value below is exact.
PANEL_POPS <- data.frame(
  pop_abb     = c("POPA", "POPB", "POPC"),
  num_subj    = c(4L, 5L, 6L),
  sup_pop_abb = c("SUP1", "SUP1", "SUP2"),
  stringsAsFactors = FALSE
)

PANEL_GENO <- list(
  # rsid,  chr, bp,   POPA,   POPB,    POPC
  list("rs1", 1L, 1000L, "0120", "01201", "012012"),
  list("rs2", 1L, 1100L, "1111", "11111", "111111"),  # no variance -> sim_z NA
  list("rs3", 1L, 1200L, "0000", "00000", "000000"),  # monomorphic -> sim_z NA
  list("rs4", 2L, 1300L, "2100", "21002", "210021"),
  list("rs5", 2L, 1400L, "0121", "01210", "012101"),
  # POPA fixed at 2 and POPC at 0, so the simulated AF1 of this SNP is exactly
  # n_POPA / (n_POPA + n_POPC): it reports which populations were actually
  # sampled, and in what proportion.
  list("rs6", 2L, 1500L, "2222", "01201", "000000")
)

# rsid of the SNP above, for tests that need to read back the sampled mix.
PANEL_DISCRIMINATING_SNP <- "rs6"

# Builds the panel in `dir` and returns the paths plus the expected values that
# the tests assert against.
make_test_panel <- function(dir = tempfile()) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  geno_file <- file.path(dir, "geno.gz")
  idx_file  <- file.path(dir, "index.gz")
  pop_file  <- file.path(dir, "pop_desc.txt")

  # Genotype line: one string per population, then that population's AF1.
  af1_pop <- function(s) sum(as.integer(strsplit(s, "")[[1]])) / (2 * nchar(s))
  geno_lines <- vapply(PANEL_GENO, function(snp) {
    g <- unlist(snp[4:6])
    paste0(paste(g, collapse = " "), " ",
           paste(sprintf("%.5f", vapply(g, af1_pop, numeric(1))), collapse = " "),
           "\n")
  }, character(1))

  fpos <- write_bgzf(geno_file, geno_lines)

  # Overall AF1 across all populations combined, as cal_af1ref() computes it.
  af1ref <- vapply(PANEL_GENO, function(snp) {
    g <- paste0(unlist(snp[4:6]), collapse = "")
    sum(as.integer(strsplit(g, "")[[1]])) / (2 * nchar(g))
  }, numeric(1))

  # Index line: rsid chr bp a1 a2 af1ref fpos
  index_lines <- vapply(seq_along(PANEL_GENO), function(i) {
    snp <- PANEL_GENO[[i]]
    sprintf("%s %d %d A G %.5f %s\n", snp[[1]], snp[[2]], snp[[3]],
            round(af1ref[i], 5), format(fpos[i], scientific = FALSE))
  }, character(1))
  write_bgzf(idx_file, index_lines, lines_per_block = 3L)

  write.table(PANEL_POPS, pop_file, quote = FALSE, sep = " ",
              row.names = FALSE, col.names = TRUE)

  list(
    dir        = dir,
    geno       = geno_file,
    index      = idx_file,
    pop_desc   = pop_file,
    geno_lines = sub("\n$", "", geno_lines),
    af1ref     = af1ref,
    rsid       = vapply(PANEL_GENO, function(s) s[[1]], character(1)),
    chr        = vapply(PANEL_GENO, function(s) s[[2]], integer(1)),
    bp         = vapply(PANEL_GENO, function(s) s[[3]], integer(1)),
    num_pops   = nrow(PANEL_POPS)
  )
}

# A description file that disagrees with the genotype data, for the bounds check.
write_bad_pop_desc <- function(dir, num_subj = c(500L, 5L, 6L)) {
  path <- file.path(dir, "bad_pop_desc.txt")
  d <- PANEL_POPS; d$num_subj <- num_subj
  write.table(d, path, quote = FALSE, sep = " ",
              row.names = FALSE, col.names = TRUE)
  path
}
