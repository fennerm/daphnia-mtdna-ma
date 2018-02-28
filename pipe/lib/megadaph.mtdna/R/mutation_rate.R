#' Calculate the mutation rate from variant and line info tables
#' @importFrom dplyr "%>%" summarize
#' @export
calc_mutation_rate <- function(variant_table, line_info) {
  mean_gen <- line_info %>% summarize(mean(generations))
  sum_bp <- line_info %>% summarize(sum(bp))
  mu <- mutation_rate(variant_table$af, mean_gen, sum_bp)
  mu <- unlist(mu)
  names(mu) <- NULL
  mu
}

#' Bootstrap difference in mean allele frequency between two groups
#'
#' See:
#' https://stats.stackexchange.com/questions/136661
#' Also:
#' Efron's and Tibshirani's - intro to the bootstrap, page 223
#' @param variant_table A variant table
#' @param line_info A line info table
#' @param group1,group2 Character; Group IDs in variant_table
#' @param by Column name in variant_table in which groups can be found
#' @param reps integer; Number of bootstrap replicates
#' @importFrom fen.R.util select_groups
#' @importFrom dplyr filter summarize
#' @return A p-value
#' @export
bootstrap_allele_frequency_test <- function(
  variant_table, line_info, by, group1, group2, reps = 1e5) {
  groups <- c(group1, group2)
  line_info_subset <- select_groups(line_info, by, groups)
  variant_table_subset <- select_groups(variant_table, by, groups)

  # Calculate the overall mean allele frequency across both groups
  combined_mean_af <- mean(variant_table_subset$af)
  observed_mean_afs <- calc_mean_af_by_group(
    variant_table = variant_table_subset, line_info = line_info_subset, by = by)

  # Construct bootstrapped mutation rate sample for each group
  bootstrapped_null_dists <- lapply(groups, function(group) {
      variant_table_subset <- select_groups(variant_table, by, group)
      line_info_subset <- select_groups(line_info, by, group)
      bootstrap_centered_null(
        variant_table = variant_table_subset,
        line_info = line_info_subset,
        statistic = indexed_mean,
        center = combined_mu,
        reps = reps)
  })

  # Calculate p-value
  p <- calculate_bootstrap_diff_p(observed_mean_afs, bootstrapped_null_dists)
  p
}

#' Bootstrap difference in mutation rate between two groups
#'
#' See:
#' https://stats.stackexchange.com/questions/136661
#' Also:
#' Efron's and Tibshirani's - intro to the bootstrap, page 223
#' @param variant_table A variant table
#' @param line_info A line info table
#' @param group1,group2 Character; Group IDs in variant_table
#' @param by Column name in variant_table in which groups can be found
#' @param reps integer; Number of bootstrap replicates
#' @importFrom fen.R.util filter select_groups
#' @importFrom dplyr summarize
#' @return A p-value
#' @export
bootstrap_mutation_rate_test <- function(
  variant_table, line_info, by, group1, group2, reps = 1e5) {
  groups <- c(group1, group2)
  line_info <- select_groups(line_info, by, groups)
  variant_table <- select_groups(variant_table, by, groups)

  # Calculate the overall mutation rate across both groups
  combined_mu <- calc_mutation_rate(variant_table_subset, line_info_subset)

  observed_mu <- mapply(
    calc_mutation_rate,
    split_table(variant_table, by),
    split_table(line_info, by))

  # Construct bootstrapped mutation rate sample for each group
  bootstrapped_null_dists <- lapply(groups, function(group) {
      variant_table_subset <- select_groups(variant_table, by, group)
      line_info_subset <- select_groups(line_info, by, group)
      bootstrap_centered_null(
        variant_table = variant_table_subset,
        line_info = line_info_subset,
        statistic = indexed_mutation_rate,
        center = combined_mu,
        gen = line_info_subset$generations,
        nuc = line_info_subset$bp,
        reps = reps)
  })
  p <- calculate_bootstrap_diff_p(observed_mu, bootstrapped_null_dists)
  p
}

#' Apply a two group statistical test to each group in the tables
#' @importFrom fen.R.util split_table ulapply
#' @export
boot_compare_all <- function(
  variant_table, line_info, test, by, within = NULL, reps = 1e5) {

  ## Bootstrap compare across groups in 'by'
  compare_within_groups <- function(variant_table, line_info) {
    by_levels <- unique(line_info[, by])

    pairwise_comparisons <- combn(by_levels, 2)
    ps <- apply(pairwise_comparisons, 2, function(x) {
      test(
        variant_table = variant_table,
        line_info = line_info,
        by = by,
        group1 = x[1],
        group2 = x[2],
        reps = reps)})
    names <- apply(pairwise_comparisons, 2, paste0, collapse = "x")
    names(ps) <- names
    ps
  }

  # Split the table by 'within'
  line_info_tables <- split_table(line_info, by = within)
  variant_tables <- split_table(variant_table, by = within)

  # Apply the test to each group in within separately
  comparisons <- unlist(mapply(
    compare_within_groups,
    variant_tables,
    line_info_tables))
  comparisons
}

#' Bootstrap the mutation rate with confidence intervals
#' @importFrom fen.R.util split_table ulapply
boot_mu_with_quantiles <- function(
  variant_table, line_info, by = NULL, reps = 1e5) {

  bootstrap_single_group <- function(variant_table_subset, line_info_subset) {
    mu <- calc_mutation_rate(
      variant_table = variant_table_subset,
      line_info = line_info_subset)
    quantiles <- boot_quantiles(
      variant_table = variant_table_subset,
      line_info = line_info_subset,
      statistic = indexed_mutation_rate,
      gen = line_info_subset$generations,
      nuc = line_info_subset$bp,
      reps = reps)
    c(mu, quantiles)
  }

  tables <- lapply(list(variant_table, line_info), split_table, by = by)

  mu_with_ci <- mapply(
    bootstrap_single_group,
    variant_table = tables[[1]],
    line_info = tables[[2]])
  species <- ulapply(tables[[2]], function(x) unique(x$species))
  species <- as.data.frame(species, stringsAsFactors = FALSE)
  groups <- colnames(mu_with_ci)
  mu_with_ci <- cbind(species, groups, t(mu_with_ci))
  colnames(mu_with_ci) <- c("species", "group", "value", "ci1", "q25", "q75", 
                            "ci2")
  rownames(mu_with_ci) <- NULL
  mu_with_ci
}


# ==============================================================================
# Helper functions
# ==============================================================================

#' Get mutant_allele frequencies for a sample
#' @import dplyr
get_afs_by_sample <- function(x, variant_table) {
  unlist(filter(variant_table, sample == x) %>% select(af))
}

#' Center a distribution at a given mean value
#' @param d numeric; the distribution
#' @param target numeric; the new mean
#' @return The centered distribution
center_distribution <- function(d, target) {
  d - mean(d) + target
}

#' Bootstrap the null distribution for a statistic centered at a given value
#'
#' This function generates a distribution for the null hypothesis that two
#' groups do not differ in some statistic. For example if we were interested in
#' whether two samples had a different mean; this function would produce a
#' bootstrapped sample for one group centered at the actual mean across both
#' groups. This ensures that the bootstrapped distribution corresponds to the
#' actual null hypothesis. Explained better here:
#' https://stats.stackexchange.com/questions/136661
#' @importFrom fen.R.util select_groups
#' @importFrom boot boot
bootstrap_centered_null <- function(
  variant_table, line_info, statistic, center, ..., reps = 1e5) {

  # Get mutation allele frequencies split by sample
  af_by_sample <- lapply(
    line_info$sample,
    get_afs_by_sample,
    variant_table = variant_table)

  # Bootstrap the mutation rate under the null
  boot_sample <- boot(
    data = af_by_sample,
    statistic = statistic,
    R = reps, ...)

  # Center the bootstrap sample at the overall mutation rate to satisfy
  # the null hypothesis
  centered_boot <- center_distribution(boot_sample$t[, 1], center)
  centered_boot
}

#' Calculate the p-value for the bootstrap difference test
#' @param observations Numeric; Vector of length == 2 giving the observed value
#'                     of the statistic on each group
#' @param bootstrapped_null list; List of length == 2 giving the centered
#'                          bootstrapped null distributions for each group
#' @return Numeric; A p-value
calculate_bootstrap_diff_p <- function(observations, bootstrapped_null) {
  observed_diff <- abs(observations[1] - observations[2])

  bootstrapped_diff <- abs(bootstrapped_null[[1]] - bootstrapped_null[[2]])

  p <- (length(which(bootstrapped_diff >= observed_diff)) + 1) /
    (length(bootstrapped_diff))
  p
}

#' Calculate the haploid mutation rate
#' @param muts numeric; mutant allele frequencies
#' @param gen numeric; Mean generation number
#' @param nuc integer; Number of nucleotides surveyed
#' @return The mutation rate
mutation_rate <- function(muts, gen, nuc) {
  arg_lengths <- c(length(gen), length(nuc))
  if (any(arg_lengths > 1)) {
    stop("Invalid gen or nuc value")
  } else if (any(c(gen, nuc) == 0)) {
    0
  } else {
    unlist(sum(muts) / (gen * nuc))
  }
}

#' Calculate the mean allele frequency for each group in by
#' @importFrom fen.R.util select_groups ulapply
calc_mean_af_by_group <- function(variant_table, line_info, by) {
  levels <- unique(line_info[, by])
  mean_afs <- ulapply(levels, function(level) {
    variant_table_subset <- select_groups(variant_table, by, level)
    mean(variant_table_subset$af)})
  mean_afs
}

#' Bootstrap statistic; take the mean
indexed_mean <- function(data, indices) {
  sample_data <- unlist(data[indices])
  c(mean(sample_data), var(sample_data))
}

#' Bootstrap statistic; calculate the mutation rate
indexed_mutation_rate <- function(data, indices, gen, nuc) {
  mutation_rate <- mutation_rate(
    unlist(data[indices]), mean(gen[indices]), sum(nuc[indices]))
  mutation_rate_by_sample <- lapply(indices, function(i) {
    mutation_rate(data[[i]], gen[i], nuc[i])})
  mutation_rate_variance <- var(unlist(mutation_rate_by_sample))
  c(mutation_rate, mutation_rate_variance)
}

#' Calculate a given statistic and compute bias corrected confidence intervals
#' @importFrom boot boot boot.ci
#' @importFrom fen.R.util select_groups
boot_quantiles <- function(variant_table, line_info, statistic, ...,
                                     reps = 1e5) {
  af_by_sample <- lapply(
    line_info$sample, get_afs_by_sample, variant_table = variant_table)
  bootstrapped_sample <- boot(
    data = af_by_sample, statistic = statistic, ..., R = reps)
  quantiles <- quantile(bootstrapped_sample$t[,1], c(0.25, 0.75))
  ci <- boot.ci(bootstrapped_sample, type = "bca")$bca[4:5]
  c(ci[1], quantiles, ci[2])
}

