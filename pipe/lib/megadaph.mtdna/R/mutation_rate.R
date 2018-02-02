variant_table <- read.csv("~/fmacrae/daphnia-mtdna-ma.private/data/mutations/magna.mt.annot.csv", stringsAsFactors = FALSE)
line_info <- read.csv("~/fmacrae/daphnia-mtdna-ma.private/daphnia-mtdna-ma/pipe/input/metadata/line_info.csv", stringsAsFactors = FALSE)
by <- "isolate"

#' Calculate the haploid mutation rate
#' @param muts numeric; mutant allele frequencies
#' @param gen numeric; Mean generation number
#' @param nuc integer; Number of nucleotides surveyed
#' @return The mutation rate
#' @export
mutation_rate <- function(muts, gen, nuc) {
  arg_lengths <- c(length(gen), length(nuc))
  if (any(arg_lengths > 1)) {
    stop("Invalid gen or nuc value")
  } else if (any(c(gen, nuc) == 0)) {
    0
  } else {
    sum(muts) / (gen * nuc)
  }
}

#' Bootstrap difference in mutation rate between two groups
#'
#' #' https://stats.stackexchange.com/questions/136661/using-bootstrap-under-h0\
#' -to-perform-a-test-for-the-difference-of-two-means-repl
#' Efron's and Tibshirani's - intro to the bootstrap, page 223
#' @param variant_table A variant table
#' @param line_info A line info table
#' @param group1,group2 Character; Group IDs in variant_table
#' @param by Column name in variant_table in which groups can be found
#' @param reps integer; Number of bootstrap replicates
#' @importFrom fen.R.util select_groups
#' @importFrom dplyr summarize
#' @return A p-value
#' @export
bootstrap_mutation_rate_test <- function(
  variant_table, line_info, by, group1, group2, reps = 1e5) {
  groups <- c(group1, group2)
  names(groups) <- groups

  # Exclude rows which are not in groups
  line_info_subset <- select_groups(line_info, by, groups)
  variant_table_subset <- select_groups(variant_table, by, groups)

  # Calculate the overall mutation rate across both groups
  combined_mu <- calc_mutation_rate_from_tables(
    variant_table_subset, line_info_subset)

  # Construct bootstrapped mutation rate sample for each group
  bootstrapped_groups <- lapply(
    groups, bootstrap_mutation_rate_by, by = by,
    variant_table = variant_table_subset, line_info = line_info_subset,
    combined_mu = combined_mu, reps = reps)

  # Observed difference in mutation rates
  observed_diff <- abs(bootstrapped_groups[[1]]$mu -
    bootstrapped_groups[[2]]$mu)

  bootstrapped_diff <- abs(bootstrapped_groups[[1]]$boot_mu -
    bootstrapped_groups[[2]]$boot_mu)

  # Number of boot samples more extreme than observed / number of boot samples
  p <- ((length(which(bootstrapped_diff >= observed_diff))) + 1) /
    (length(bootstrapped_diff))
  p
}

#' Bootstrap mutation rate difference for all possible group comparisons in 'by'
#' @param variant_table A variant table
#' @param line_info A line info table
#' @param within Column name in variant_table; If specified the comparisons in
#'               'by' are made separately for each group in 'within'.
#' @param group1,group2 Character; Group IDs in variant_table
#' @param by Column name in variant_table in which groups can be found
#' @importFrom fen.R.util select_groups
#' @export
compare_all_mutation_rates <- function(
  variant_table, line_info, by, within = NULL, reps = 1e5) {
  if (!is.null(within)) {
    # Split the inputs by within then call helper function on each seperately.
    within_levels <- unique(line_info[, within])

    bootstrapped_p <- lapply(within_levels, function(within_level) {
      variant_table_within <- select_groups(variant_table, within, within_level)
      line_info_within <- select_groups(line_info, within, within_level)
      compare_mutation_rates_within(
        variant_table = variant_table_within, line_info = line_info_within,
        by = by, reps = reps)})
    bootstrapped_p <- unlist(bootstrapped_p, recursive = FALSE)

  } else {
    # Just call the helper function on the entire input
    bootstrapped_p <- compare_mutation_rates_within(
      variant_table = variant_table_within, line_info = line_info_within,
      by = by, reps = reps)
  }
  bootstrapped_p
}

#' Bootstrap the mutation rate across samples and return the quantiles
#' @param af list of numeric; mutant allele frequencies split by sample
#' @param gen list of numeric; mean generation number split by sample
#' @param nuc list of integer; bases surveyed split by sample
#' @return Numeric vector; First item is the mutation rate. The
#'         other five are the 2.5%, 25%, 50%, 75%, and 97.5% quantiles.
#' @export
bootstrap_mutation_rate_quantiles <- function(af, gen, nuc, reps=1e5) {
  mu <- mutation_rate(unlist(af), mean(gen), sum(nuc))
  boot <- bootstrap_mutation_rate(af, gen, nuc, reps)
  quantiles <- quantile(boot, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  c(mu, quantiles)
}

#' @importFrom fen.R.util select_groups
bootstrap_mutation_rate_quantiles_from_tables <- function(
  variant_table, line_info, by = NULL) {
  by_levels <- unique(line_info[, by])
  bootstrapped_quantiles <- sapply(by_levels, function(lev) {
    variant_table_subset <- select_groups(variant_table, by, lev)
    line_info_subset <- select_groups(line_info, by, lev)
    af <- lapply(line_info_subset[, "sample"], get_afs_by_sample,
      variant_table_subset)
    bootstrap_mutation_rate_quantiles_from_tables(af, gen, nuc, reps)
  })
}

#' Bootstrap the haploid mutation rate statistic across samples
#' @param af list of numeric; mutant allele frequencies split by sample
#' @param gen list of numeric; mean generation number split by sample
#' @param nuc list of integer; bases surveyed split by sample
#' @param reps integer; Number of bootstrap replicates
#' @return A bootstrapped sample of mutation rates
#' @export
bootstrap_mutation_rate <- function(af, gen, nuc, reps = 1e5) {
  nsamples <- length(gen)
  sample_ids <- 1:nsamples
  sum_nuc <- sum(nuc)
  bt <- replicate(reps, {
    boot_sample <- sample(sample_ids, replace = TRUE)
    mutation_rate(unlist(af[boot_sample]), mean(gen[boot_sample]), sum_nuc)},
    simplify = FALSE)
  unlist(bt)
}

#' Calculate the mutation rate from variant and line info tables
#' @import dplyr
#' @export
calc_mutation_rate_from_tables <- function(variant_table, line_info) {
  mean_gen <- line_info %>% summarize(mean(generations))
  sum_bp <- line_info %>% summarize(sum(bp))
  mu <- mutation_rate(variant_table$af, mean_gen, sum_bp)
  unlist(mu)
}

# ==============================================================================
# Helper functions
# ==============================================================================


#' Get mutant_allele frequencies for a sample
#' @import dplyr
get_afs_by_sample <- function(sample, variant_table) {
  unlist(filter(variant_table, sample == sample) %>% select(af))
}

#' Center a distribution at a given mean value
#' @param d numeric; the distribution
#' @param target numeric; the new mean
#' @return The centered distribution
center_distribution <- function(d, target) {
  d - mean(d) + target
}

#' Bootstrap the mutation rate for a specific group
#' @importFrom fen.R.util select_groups
bootstrap_mutation_rate_by <- function(
  group, by, variant_table, line_info, combined_mu, reps) {
  # Extract group from variant_table and line_info_subset
  group_variant_table <- select_groups(variant_table, by, group)
  group_line_info <- select_groups(line_info, by, group)

  # Get mutation allele frequencies split by sample
  mutant_afs <- lapply(group_line_info$sample, get_afs_by_sample,
    variant_table = group_variant_table)

  mu <- calc_mutation_rate_from_tables(group_variant_table, group_line_info)

  # Bootstrap the mutation rate under the null
  boot_mu <- bootstrap_mutation_rate(
    mutant_afs, group_line_info$generations, group_line_info$bp, reps)

  # Center the bootstrap sample at the overall mutation rate to satisfy
  # the null hypothesis
  centered_boot_mu <- center_distribution(boot_mu, combined_mu)
  list(mu = mu, boot_mu = centered_boot_mu)
}

#' Helper function called by 'compare_all_mutation_rates'. This function in turn
#' just repeatedly calls the 2-group bootstrap test function for each pairwise
#' combination of levels in 'by'
compare_mutation_rates_within <- function(
  variant_table, line_info, by, reps) {
  by_levels <- unique(line_info[, by])

  # Call bootstrap test on each possible pairwise combination of the grouping
  pairwise_comparisons <- combn(by_levels, 2)
  d <- apply(pairwise_comparisons, 2, function(x) {
    bootstrap_mutation_rate_test(
      variant_table = variant_table, line_info = line_info, by = by,
      group1 = x[1], group2 = x[2], reps = reps)[1]})
  names <- apply(pairwise_comparisons, 2, paste0, collapse = "")
  names(d) <- names
  d
}
