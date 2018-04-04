#' Bootstrap statistic; take the mean of the mutant allele frequencies
#' @export
indexed_mean_af <- function(data, indices) {
  afs <- unlist(data[indices, "af"])
  c(mean(afs), var(afs))
}

#' Bootstrap statistic; calculate the mutation rate
#' @importFrom magrittr set_names "%>%"
#' @importFrom dplyr slice summarize ungroup
#' @export
indexed_mutation_rate <- function(dat, indices) {
  dat_subset <- dat %>% ungroup %>% slice(indices)
  mu <- dat_subset %>% summarize(mutation_rate(af, mean(generations), sum(bp)))
  mu_by_sample <- dat_subset %>%
    rowwise %>%
    summarize(mutation_rate(af, generations, bp))
  var <- var(unlist(mu_by_sample))
  result <- unlist(c(mu, var)) %>% set_names(c('mutation_rate', 'variance'))
  result
}


#' Calculate the haploid mutation rate
#' @param muts numeric; mutant allele frequencies per sample
#' @param gen numeric; generation number per sample
#' @param bp integer; Number of nucleotides surveyed per sample
#' @return The mutation rate
mutation_rate <- function(muts, gen, bp) {
  muts <- unlist(muts)
  gen <- unlist(gen)
  bp <- unlist(bp)
  sum(muts) / (mean(gen) * sum(bp))
}
