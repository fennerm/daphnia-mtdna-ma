# Creating and manipulating nested variant tables.

#' Merge table of mutations with the metadata table.
#'
#' @param table_of_mutations data.frame; A table of discovered mutations,
#'   containing at least the following columns: 'sample', 'pos', 'ref', 'alt',
#'   'class' (one of 'snv', 'insertion', 'deletion'). If mutations are
#'   mitochondrial the table should also contain an 'af' column with mutant
#'   allele frequencies.
#' @param line_info data.frame; A table of group metadata. E.g 'population',
#'   'genotype' etc. One column must contain samples (same as in
#'   `mutation_table`. All info from this table will be added to the mutation
#'   table.
#' @return A nested tibble
#' @importFrom magrittr "%>%"
#' @importFrom dplyr funs group_by right_join summarize_all vars
#' @importFrom tidyr replace_na nest
#' @export
VariantTable <- function(table_of_mutations, line_info) {
  names(table_of_mutations) <- tolower(names(table_of_mutations))
  list_cnames <- names(table_of_mutations)[names(table_of_mutations) != "sample"]
  nonlist_cnames <- names(line_info)[names(line_info) != "sample"]
  variant_table <- table_of_mutations %>%
      group_by(sample) %>%
      # Merge table_of_mutations with line_info
      right_join(., line_info, by = "sample") %>%
      # Convert certain columns to list
      mutate_at(vars(list_cnames), list) %>%
      # Converting columns to lists creates duplicate rows, remove them
      nest(.key = "tmp") %>%
      mutate(data = tmp %>% map(~.[1,])) %>%
      select(sample, data)
   variant_table
}

#' Return groupings which vary between samples but not within
#'
#' @importFrom magrittr "%>%"
#' @importFrom dplyr summarize_all
get_per_sample_levels <- function(variant_table) {
  colclasses <- variant_table %>% summarize_all(class)
  is_per_sample_level <- colclasses != "list" & names(variant_table) != "sample"
  names(variant_table)[which(is_per_sample_level)]
}


#' Change the grouping of a variant table.
#'
#' @param variant_table tibble Returned by `VariantTable` constructor
#' @param by character Main grouping column. If by == "combined", the table will
#'   nested as a single group.
#' @param within character Inner grouping column
#' @return tibble
#' @export
regroup <- function(variant_table, by) {
  if (by == "combined") {
    variant_table$combined = "combined"
  }
  regrouped_table <- variant_table %>%
    unnest %>%
    group_by(!!as.name(by)) %>%
    nest
  regrouped_table
}

#' Determine if variants in table are from mitochondrial or nuclear genome
#'
#' @return One of "mt" or "nuc"
get_genome <- function(variant_table) {
  column_names <- names(variant_table)
  if ('af' %in% column_names) {
    return("mt")
  } else if ('data' %in% column_names) {
    if ('af' %in% names(variant_table$data[[1]])) {
      return("mt")
    }
  }
  return("nuc")
}
