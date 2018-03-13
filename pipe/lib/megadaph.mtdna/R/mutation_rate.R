#' Bootstrap statistic; take the mean of the mutant allele frequencies
#' @export
indexed_mean_af <- function(data, indices) {
  afs <- unlist(data[indices, "afs"])
  c(mean(afs), var(afs))
}

#' Bootstrap statistic; calculate the mutation rate
#' @importFrom magrittr set_names "%>%"
#' @importFrom dplyr slice summarize ungroup
#' @export
indexed_mutation_rate <- function(dat, indices) {
  dat_subset <- dat %>% ungroup %>% slice(indices)
  mu <- dat_subset %>% summarize(mutation_rate(afs, generations, bp))
  mu_by_sample <- dat_subset %>%
    rowwise %>%
    summarize(mutation_rate(afs, generations, bp))
  var <- var(unlist(mu_by_sample))
  result <- unlist(c(mu, var)) %>% set_names(c('mutation_rate', 'variance'))
  result
}


#' Get mutant_allele frequencies for a sample
#' @import dplyr
get_afs_by_sample <- function(x, variant_table) {
  afs <- unlist(filter(variant_table, sample == x) %>% select(af))
  if (length(afs) == 0) {
    afs <- 0
  }
  afs
}


#' Calculate the haploid mutation rate
#' @param muts numeric; mutant allele frequencies
#' @param gen numeric; Mean generation number
#' @param bp integer; Number of nucleotides surveyed
#' @return The mutation rate
mutation_rate <- function(muts, gen, bp) {
  muts <- unlist(muts)
  gen <- unlist(gen)
  bp <- unlist(bp)
  arg_lengths <- c(length(gen), length(bp))
  if (any(arg_lengths > 1)) {
    stop("Invalid gen or bp value: must by single number")
  } else if (any(c(gen, bp) == 0)) {
    0
  } else {
    sum(muts) / (gen * bp)
  }
}
