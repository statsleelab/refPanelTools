#' Extract Genomic Region Data as a Data Frame
#'
#' Returns the SNPs in a base-pair range as a data.frame, combining metadata
#' from the index file with genotype strings and per-population allele
#' frequencies.
#'
#' Metadata and genotype records are collected together in a single pass over
#' the index, so the cost is proportional to the size of the region rather than
#' the size of the index.
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
#'     \item{geno_<pop>}{Character. Genotype string for each population (0/1/2
#'       per subject, packed as a single string; leading zeros are preserved).}
#'     \item{af1_<pop>}{Numeric. Alternate allele frequency for each population.}
#'   }
#'   A region holding no SNPs gives a data.frame with zero rows and these same
#'   columns, so results from several regions can be combined with `rbind()`
#'   without special-casing the empty result.
#'
#' @seealso [extract_reg_data()] to write the genotype records to a file
#'   instead, which is the better choice for a region too large to hold in
#'   memory.
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

  rec <- read_reg_records(chr_num, start_bp, end_bp, num_pops,
                          index_data_file, reference_data_file)

  geno <- as.data.frame(rec$geno, stringsAsFactors = FALSE)
  af1  <- as.data.frame(rec$af1)
  names(geno) <- paste0("geno_", pop_names)
  names(af1)  <- paste0("af1_",  pop_names)

  cbind(
    data.frame(rsid   = rec$rsid,
               chr    = rec$chr,
               bp     = rec$bp,
               a1     = rec$a1,
               a2     = rec$a2,
               af1ref = rec$af1ref,
               stringsAsFactors = FALSE),
    geno, af1
  )
}
