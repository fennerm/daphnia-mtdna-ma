#' Add q-values to a list of test tables
#' @param test_tables List of test tables
#' @param fdr_level False discovery rate to control by
#' @return List of test tables
#' @export
#' @importFrom qvalue qvalue
mult_comparisons_correct <- function(test_tables, fdr_level=0.005) {
  for (i in c(1:length(test_tables))) {
    test_tables[[i]][, "mult_comp_group"] <- i
  }
  # Merge into a single table
  merged_table <- do.call("rbind", test_tables)

  # Calculate q values
  fdr <- qvalue(merged_table$p_value, fdr.level = fdr_level)

  # Add the FDR information to the table
  fdr_table <- cbind(merged_table, qvalue = fdr$qvalues,
                     significant = fdr$significant)
  colnames(fdr_table)
  fdr_table <- data.frame(fdr_table, stringsAsFactors = FALSE)

  # Split the table back into genotypes
  fdr_tables <- split(fdr_table, fdr_table$mult_comp_group)
  for (i in c(1:length(fdr_tables))) {
    fdr_tables[[i]]$mult_comp_group <- NULL
  }
  fdr_tables
}
