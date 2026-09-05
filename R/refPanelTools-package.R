#' @keywords internal
#'
#' @details
#' The package indexes BGZF-compressed reference panels, extracts genotype data
#' by chromosome, genomic region or population subset, and simulates allele
#' frequencies and association Z-scores under the null for multi-ethnic GWAS.
#'
#' @seealso
#' [indexer()], [extract_reg_data()], [read_region()], [simulate_af1_z()]
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib refPanelTools, .registration = TRUE
"_PACKAGE"
