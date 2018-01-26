#' Exclude nonsignificant, non-unique and high-strand-bias variants
#' @param fdr_table A test table with qvalues
#' @param max_strand_bias Numeric; Maximum strand bias allowed for included
#'                        variants
#' @return A filtered variant table
#' @importFrom dplyr filter
#'@export
filter_variants <- function(fdr_table, max_strand_bias) {
  filter(fdr_table, fdr_table$unique, fdr_table$significant,
         fdr_table$strand_bias < max_strand_bias)
}
