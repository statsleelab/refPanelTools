#' Extract Genomic Region Data as a Data Frame
#'
#' A convenience wrapper around [extract_reg_data()] that returns the extracted
#' data as a data.frame, combining SNP metadata from the index file with
#' genotype strings and per-population allele frequencies.
#'
#' The raw output of [extract_reg_data()] contains one row per SNP with
#' space-separated genotype strings (one per population) followed by
#' per-population AF1 values. This wrapper parses that output and attaches
#' SNP metadata (rsid, chr, bp, a1, a2) from the index file.
#'
#' @param chr_num Integer. Chromosome number (1--22).
#' @param start_bp Integer. Start base pair position (inclusive).
#' @param end_bp Integer. End base pair position (inclusive).
#' @param index_data_file Character. Path to the BGZF-compressed reference
#'   panel index file, whose columns are `rsid chr bp a1 a2 af1ref fpos`.
#' @param reference_data_file Character. Path to the BGZF-compressed reference
#'   genotype file.
#' @param pop_desc_file Character. Path to the population description file
#'   (columns: `pop_abb num_subj sup_pop_abb`, with a header row). Used to
#'   label genotype and AF1 columns. If `NULL`, `num_pops` must be supplied
#'   and generic names (`geno_pop1`, `af1_pop1`, ...) are used.
#' @param num_pops Integer. Total number of populations in the reference panel.
#'   Required when `pop_desc_file` is `NULL`; ignored otherwise.
#'
#' @return A data.frame with one row per SNP in the region and columns:
#'   \describe{
#'     \item{rsid}{SNP identifier.}
#'     \item{chr}{Chromosome number.}
#'     \item{bp}{Base pair position.}
#'     \item{a1, a2}{Reference and alternate alleles.}
#'     \item{af1ref}{Overall reference panel AF1.}
#'     \item{geno_<pop>}{Genotype string for each population (0/1/2 per subject,
#'       packed as a single character string).}
#'     \item{af1_<pop>}{Allele frequency of alternate allele for each population.}
#'   }
#'   Returns an empty data.frame if no SNPs fall in the requested region.
#'
#' @seealso [extract_reg_data()] for the underlying file-based function.
#'
#' @importFrom utils read.table
#' @examples
#' \dontrun{
#' df <- read_region(
#'   chr_num             = 14,
#'   start_bp            = 104000000,
#'   end_bp              = 104200000,
#'   index_data_file     = "/33KG/33kg_index.gz",
#'   reference_data_file = "/33KG/33kg_geno.gz",
#'   pop_desc_file       = "/33KG/33kg_pop_desc.txt"
#' )
#' head(df[, c("rsid", "chr", "bp", "a1", "a2", "af1ref")])
#' }
#' @export
read_region <- function(chr_num, start_bp, end_bp,
                        index_data_file, reference_data_file,
                        pop_desc_file = NULL, num_pops = NULL) {

  # Resolve population names and count
  if (!is.null(pop_desc_file)) {
    pop_desc  <- read.table(pop_desc_file, header = TRUE,
                            stringsAsFactors = FALSE)
    pop_names <- as.character(pop_desc[[1]])
    num_pops  <- length(pop_names)
  } else {
    if (is.null(num_pops))
      stop("Provide either 'pop_desc_file' or 'num_pops'.")
    pop_names <- paste0("pop", seq_len(num_pops))
  }

  # Extract genotype data to a temp file, then read it back
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)

  extract_reg_data(chr_num, start_bp, end_bp, num_pops,
                   index_data_file, reference_data_file, tmp)

  if (file.info(tmp)$size == 0L)
    return(data.frame())

  geno_data <- read.table(tmp, header = FALSE, sep = "",
                          colClasses = "character")

  colnames(geno_data) <- c(paste0("geno_", pop_names),
                            paste0("af1_",  pop_names))

  # Read SNP metadata from the index file and filter to the requested region
  # Tip: data.table::fread() is significantly faster for large index files
  idx <- read.table(index_data_file, header = FALSE, sep = "",
                    col.names   = c("rsid", "chr", "bp",
                                    "a1", "a2", "af1ref", "fpos"),
                    colClasses  = c("character", "integer", "integer",
                                    "character", "character", "numeric",
                                    "numeric"))
  snp_meta <- idx[idx$chr == chr_num &
                  idx$bp  >= start_bp &
                  idx$bp  <= end_bp,
                  c("rsid", "chr", "bp", "a1", "a2", "af1ref"),
                  drop = FALSE]
  rownames(snp_meta) <- NULL

  if (nrow(snp_meta) != nrow(geno_data))
    stop("Row count mismatch between index (", nrow(snp_meta),
         ") and genotype data (", nrow(geno_data),
         "). The index and genotype files may be out of sync.")

  cbind(snp_meta, geno_data)
}
