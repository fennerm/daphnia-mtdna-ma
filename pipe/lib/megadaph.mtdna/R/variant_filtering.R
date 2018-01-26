#' Exclude nonsignificant, non-unique and high-strand-bias variants
#' @param fdr_table A test table with qvalues
#' @param max_strand_bias Numeric; Maximum strand bias allowed for included
#'                        variants
#' @return A filtered variant table
#'@export
filter_variants <- function(fdr_table, max_strand_bias) {
  var_table <- fdr_table[which(fdr_table$unique), ]
  var_table <- var_table[which(var_table$significant), ]
  var_table <- var_table[which(var_table$strand_bias < max_strand_bias), ]
  var_table
}
