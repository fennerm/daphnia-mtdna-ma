## Mutation Rate Calculations

#' Calculate the haploid mutation rate
#'
#' @param af numeric or list of numeric; mutant allele frequencies per
#'   sample
#' @param generations numeric or list of numeric; generation number per sample
#' @param bases_surveyed integer or list of integer; Number of nucleotides surveyed per
#'   sample
#' @return numeric
haploid_mutation_rate <- function(af, generations, bases_surveyed) {
  sum_af <- sum(remove_na(unlist(af)))
  mean_generations <- mean(unlist(generations))
  total_bases_surveyed <- sum(unlist(bases_surveyed))
  sum_af / (mean_generations * total_bases_surveyed)
}


#' Calculate the diploid mutation rate
#'
#' @param pos list or vector of mutant positions
#' @param generations numeric or list of numeric; generation number per sample
#' @param bases_surveyed integer or list of integer; Number of nucleotides
#'   surveyed per sample
#' @return numeric
diploid_mutation_rate <- function(pos, generations, bases_surveyed) {
  num_pos <- length(remove_na(unlist(pos)))
  mean_generations <- mean(unlist(generations))
  total_bases_surveyed <- sum(unlist(bases_surveyed))
  num_pos / (2 * mean_generations * total_bases_surveyed)
}


#' Bootstrap statistic: Calculate mean allele frequency from a variant table.
indexed_mean_af <- function(variant_table, indices) {
  afs <- unlist(variant_table[indices, "af"])
  c(af = mean(afs), variance = var(afs))
}



#' Faster alternative to slice for selecting rows from nested tbl
select_rows <- function(tbl, indices) {
  as.data.frame(ungroup(tbl))[indices, ]
}

#' Bootstrap statistic; calculate the mutation rate and variance.
#'
#' @importFrom magrittr set_names "%>%"
#' @importFrom dplyr ungroup select slice
#' @importFrom purrr pmap
indexed_mutation_rate <- function(variant_table, indices) {
  tbl_subset <- select_rows(variant_table, indices)

  if (get_genome(variant_table) == "mt") {
    mutation_rate <- haploid_mutation_rate
    mut_column <- "af"
  } else {
    mutation_rate <- diploid_mutation_rate
    mut_column <- "pos"
  }
  mu <- mutation_rate(
    tbl_subset[mut_column],
    tbl_subset$generations,
    tbl_subset$bases_surveyed
  )
  mu_by_sample <- pmap(
    tbl_subset[c(mut_column, "generations", "bases_surveyed")],
    mutation_rate
  )
  variance <- var(unlist(mu_by_sample))
  result <- unlist(c(mutation_rate = mu, variance = variance))
  result
}

#' Bootstrap the mutation rate across groups
#'
#' @param variant_table Returned table from `VariantTable`
#' @param reps integer Number of bootstrap replicates
#' @return tibble with bootstapped mutation rate per group, quantiles and
#'   confidence intervals
#' @importFrom bootr boot_quantiles
#' @export
boot_mut_rate <- function(variant_table, reps=1e4) {
  boot_quantiles(variant_table, indexed_mutation_rate, reps = reps)
}

#' Compare mutation rates via bootstrapping
#'
#' @param variant_table Returned table from `VariantTable`
#' @param reps integer Number of bootstrap replicates
#' @param within character Only compare within this grouping column
#' @return tibble with bootstapped mutation rate per group, quantiles and
#'   confidence intervals
#' @importFrom bootr boot_compare_all
#' @export
boot_compare_mut_rates <- function(variant_table, within, reps=1e4) {
  boot_compare_all(
    variant_table,
    within = within,
    statistic = indexed_mutation_rate,
    reps = reps
  )
}

#' Compare mean allele frequency via bootstrapping
#'
#' Only relevant for mitochondrial datasets.
#' @param variant_table Returned table from `VariantTable`
#' @param reps integer Number of bootstrap replicates
#' @param within character Only compare within this grouping column
#' @return tibble with bootstapped mean allele frequency per group, quantiles
#'   and confidence intervals
#' @importFrom bootr boot_compare_all
#' @export
boot_compare_mean_af <- function(variant_table, within, reps=1e4) {
  boot_compare_all(
    variant_table,
    within = within,
    statistic = indexed_mean_af,
    reps = reps
  )
}
