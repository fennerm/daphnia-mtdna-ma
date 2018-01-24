#' @export
#' @importFrom qvalue qvalue
mult_comparisons_correct <- function(test_tables, fdr_level=0.05) {
  # Merge into a single table
  merged_table <- do.call(rbind, test_tables)

  # Calculate q values
  fdr <- qvalue(merged_table$p_value, fdr.level = fdr_level)

  # Add the FDR information to the table
  fdr_table <- cbind(merged_table, fdr$qvalues, fdr$significant)
  fdr_table <- data.frame(fdr_table, stringsAsFactors = FALSE)

  # Split the table back into genotypes
  fdr_tables <- split(fdr_table, fdr_table$genotype)
  fdr_tables
}

filter_variants
  fdr_table <- fdr_table[which(fdr_table$unique),]
  fdr_table <- fdr_table[which(fdr_table$strand_bias < 60),]
  fdr_table
  # split_by_isolate <- split(merged_table, merged_table$isolate)
  # fdr <- lapply(split_by_isolate, function(x) {
  #     qvalue(x$p_value, fdr.level=fdr_level/length(unique(merged_table$isolate)))
  # })
  # significant <- unlist(lapply(fdr, function(x) x$significant))
  # q_value <- unlist(lapply(fdr, function(x) x$qvalues))
  # fdr_table <- cbind(merged_table, q_value, significant)
  # fdr_table <- data.frame(fdr_table, stringsAsFactors=FALSE)
  # fdr_table <- fdr_table[which(fdr_table$unique),]
  # fdr_table <- fdr_table[which(fdr_table$strand_bias < 60),]
  # fdr_table
}

get_stat_significant <- function(fdr_table) {
  sig_table <- fdr_table[which(fdr_table$significant),]
  row.names(sig_table) <- NULL
  sig_table
}

calc_coverages <- function(var_tables) {
  coverages <- sapply(var_tables, function(x) mean(x$coverage))
  coverages
}


#' Subsample a count vector
#' @param r An integer vector
#' @param cap Integer; Vector will be subsampled such that its sum = cap
#' @return An integer vector
#' @export
subsample <- function(r, cap) {
  r <- unlist(r)
  cov <- sum(r)
  if (cov > cap) {
    # Calculate the surplus and randomly subtract from the row
    surplus <- cov - cap
    prop <- r / sum(r)
    smp <- rmultinom(1, surplus, prop)
    new_row <- as.integer(r - smp)
  } else {
    new_row <- r
  }
  new_row
}

#' Subsample an allele count pileup
#' @param pile A data.table; allele counts at each genome position
#' @param cap Integer; Pileup will be subsampled such that its sum = cap
#' @return A pileup
#' @export
#' @importFrom data.table setDT setnames
subsample_pileup <- function(pile, cap) {
  subsampled_rows <- apply(pile, 1, subsample)
  subsampled_pile <- setDT(data.frame(t(subsampled_rows)))
  setnames(subsampled_pile, colnames(pile))
  subsampled_pile
}

#' Given a column index into an unstranded pileup, return the equiv. stranded
#' indices
.stranded_indices <- function(idx) {
  str_idx <- list(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 10), c(11, 12))
  str_idx <- matrix(unlist(str_idx[idx]), ncol = 2, byrow = TRUE)
  str_idx
}

  # piles <<- lapply(tabs, "[[", 1)
  # mut_wt_matrices <<- lapply(tabs, "[[", 2)
  # test_tables <<- lapply(tabs, "[[", 3)
  # merged_table <- do.call(rbind, test_tables)
  # p_value <- readRDS("../tables/corrado_p.Rds")
  # merged_table <- cbind(merged_table, p_value)
  # test_tables <- split(merged_table, merged_table$isolate)
  # formatted_test_tables <- lapply(test_tables, merge_significant_indels)
  # formatted_test_tables <<- lapply(formatted_test_tables, rename_small_deletions)
  # fdr_table <<- mult_comparisons_correct(formatted_test_tables,
  #                                        fdr_level=0.01)
  # var_table <<- get_stat_significant(fdr_table)
  # var_table
