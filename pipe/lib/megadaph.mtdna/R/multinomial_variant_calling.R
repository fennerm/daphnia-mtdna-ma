# ==============================================================================
# EXPORTS
# ==============================================================================

#' Create a consensus sequence from a set of nucleotide count pileups
#' @param piles List of matrices; Each list item contains allele counts at
#'        each genome position for a single sample.
#' @return Character vector; The multisample consensus sequence as a string
#' @importFrom fen.R.util mode
#' @export
create_consensus <- function(piles) {
  nsamples <- length(piles)

  # Determine the highest frequency allele at each site for each sample
  major_alleles <- lapply(piles, get_major_alleles)

  # Since all the same bases are expected to be represented in each read
  # count file, we set the number of base pairs using the first file.
  major_alleles <- matrix(unlist(major_alleles), ncol = nsamples)
  consensus <- apply(major_alleles, 1, mode)
  consensus
}

#' Create a sequence of the consensus alternate (mutant) alleles at each genome
#' position among multiple samples.
#'
#' The alternate allele is the non-consensus base with the highest frequency at
#' a locus within a single sample. These alternate alleles may be true mutants
#' or might just be sequencing error.
#' @param piles List of matrices; Each list item contains allele counts at
#'        each genome position for a single sample.
#' @param consensus character; The consensus sequence for the samples
#' @return character; The alternate consensus sequence
#' @importFrom fen.R.util ulapply
#' @export
create_mutant_consensus <- function(piles, consensus) {
  bp <- nrow(piles[[1]])

  # Create table of minor allele frequencies for each sample.
  mafs <- lapply(piles, get_minor_allele_frequencies, consensus = consensus)
  mutant_consensus <- ulapply(1:bp, function(i) {
                                print(i)
                                v <- ulapply(mafs, "[", i)
                                vnames <- lapply(mafs, names)
                                vnames[[which.max(v)]][i]
                       })
  mutant_consensus
}

#' Calculate mutant allele frequency at each genome position for a sample
#' @param pile matrix; allele counts at each genome position
#' @param mutant_consensus; An alternate allele consensus sequence produced by
#'        create_mutant_consensus
#' @return Numeric vector
#' @importFrom fen.R.util ulapply
#' @export
compute_mutant_allele_frequencies <- function(pile, mutant_consensus) {
  bp <- nrow(pile)
  ulapply(1:bp, function(i) {
            compute_allele_frequency(pile[i, ], mutant_consensus[i])
          })
}

#' Calculate strand bias for each genome position using fisher tests
#' @param stranded_piles matrix; Allele count pileups separated by strand
#' @param consensus character; The consensus sequence for the samples
#' @param mutant_consensus character; An alternate allele consensus sequence
#'        produced by create_mutant_consensus
#' @param mutant_samples numeric vector; which sample is a potential mutant for
#'        each genome position
#' @return Numeric vector of p-values
#' @export
compute_all_strand_bias <- function(stranded_piles, consensus, mutant_consensus,
                                    mutant_samples) {
  bp <- length(consensus)
  # Get the potential mutant rows from the pileups
  mutant_count_list <- lapply(1:bp, function(i) {
                                unlist(stranded_piles[[mutant_samples[i]]][i,])
                                    })
  strand_bias <- unlist(mapply(compute_strand_bias, mutant_count_list,
                               consensus, mutant_consensus, SIMPLIFY = FALSE))
  strand_bias
}

#' Convert an allele count pileup to a table of mutant vs. wild-type counts
#' @param pile A matrix; allele counts at each genome position
#' @param mutant_consensus; An alternate allele consensus sequence produced by
#'        create_mutant_consensus
#' @return An nx2 matrix, where columns are mutant and wild-type allele counts,
#'         and rows are genome positions.
#' @export
convert_to_mut_wt_counts <- function(pile, mutant_consensus) {
  mut_counts <- get_allele_counts(pile, mutant_consensus)
  cov <- get_coverage(pile)
  wt_counts <- cov - mut_counts
  cbind(mut_counts, wt_counts)
}

#' Get the allele counts for a set alleles
#' @param pile A matrix; allele counts at each genome position
#' @param alleles Character vector
#' @return Numeric vector
#' @importFrom fen.R.util ulapply
#' @export
get_allele_counts <- function(pile, alleles) {
  bp <- nrow(pile)
  counts <- ulapply(1:bp, function(i) {
                      pile[i, alleles[i]]
                                    })
  counts

}

#' Get the major allele at each genome position for a single sample
#' WARNING: Returns the first value if tied
#' @param pile matrix; allele counts at each genome position
#' @return Character; Vector of major alleles, each allele one of
#'         {A, C, G, T, +, -}
#' @export
get_major_alleles <- function(pile) {
  apply(pile, 1, get_highest_frequency_allele)
}

#' Get the alleles with the second highest frequency at each position
#' @param pile matrix; allele counts at each genome position
#' @param consensus The consensus sequence for the samples
#' @return Character; Vector of minor alleles, each allele one of
#'         {A, C, G, T, +, -}
#' @importFrom fen.R.util ulapply
#' @export
get_minor_alleles <- function(pile, consensus) {
  bp <- nrow(pile)
  ulapply(1:bp, function(i) {
            get_minor_allele(pile[i, ], consensus[i])
  })
}

#' Get sequencing coverage for each genome position
#' @param pile matrix; allele counts at each genome position
#' @return A numeric vector
#' @export
get_coverage <- function(pile) {
  # Function for getting coverage of a single sample.
  get_sample_coverage <- function(p) {
    apply(pile, 1, sum)
  }

  # If list of matrices, apply to each.
  if (is.null(ncol(pile))) {
    lapply(pile, get_sample_coverage)
  } else {
    get_sample_coverage(pile)
  }
}

#' Convert an allele count pileup to a table of mutant vs. coverage
#' @param pile matrix; allele counts at each genome position
#' @param mutant_consensus; An alternate allele consensus sequence produced by
#'        `create_mutant_consensus`
#' @return An nx2 matrix, where cols are mutant allele counts, and coverage, and
#' rows are genome positions.
#' @export
convert_to_mut_cov_counts <- function(pile, mutant_consensus) {
  mut_counts <- get_allele_counts(pile, mutant_consensus)
  depths <- get_coverage(pile)
  cbind(mut_counts, depths)
}

#' Determine which sites are potientially mutant in only one sample
#' @param mut_cov_counts List of matrices; Each matrix has two columns,
#'                       first column is mutant allele counts, second column is
#'                       total sequencing coverage. Rows are genome positions
#' @param seq_error_rates Numeric vector; Sequencing error rate estimates for
#'                        each sample
#' @return Logical vector
#' @export
determine_unique <- function(mut_cov_counts, seq_error_rates) {
  nsamples <- length(mut_cov_counts)

  # Run binomial test on each genome position for each sample
  p <- mapply(binom_tests, mut_cov_counts, seq_error_rates, SIMPLIFY = FALSE)

  # Convert list to multisample matrix
  p_matrix <- matrix(unlist(p), ncol = nsamples)

  unique <- apply(p_matrix, 1, is_unique)
  unique
}

#' Get the type of mutation
#' @param allele Character; Vector of alleles, {A, C, G, T, +, -}
#' @return Character; Vector of {snv, insertion, deletion}
#' @export
get_mutation_class <- function(mutant_allele) {
  ifelse(mutant_allele == "-", "deletion",
         ifelse(mutant_allele == "+", "insertion", "snv"))
}

#' Convert a list of tables separated by sample to a list of tables separated by
#' genome position
#' @param x List of matrices or vectors
#' @return List of matrices
#' @export
by_position <- function(x) {
  # If x is a list of vectors
  if (is.null(dim(x[[1]])[2])) {
    nsamples <- length(x)
    bp <- length(x[[1]])
    y <- lapply(1:bp, function(i) ulapply(1:nsamples, function(j) x[[j]][i]))
  # If x is a list of matrices
  } else {
    bp <- nrow(x[[1]])
    y <- lapply(1:bp, function(i) {
             mat <- sapply(x, function(x) x[i, ])
             colnames(mat) <- NULL
             t(mat)
           })
  }
  y
}

#' Assign a p-value to a variant
#' @param mut_cov_matrix matrix; Rows are genome samples, column 1 is
#'                       mutant allele count. Column 2 is sequencing coverage.
#' @return A p-value
#' @importFrom pmultinom pmultinom
#' @importFrom fen.R.util ulapply
#' @export
call_variant <- function(mut_cov_matrix) {
  # Model mutation as a multinomial distribution. Under the null hypothesis, the
  # chance of observing a mutant allele is determined by sequencing coverage
  # only.

  mut_af <- mut_cov_matrix[, 1] / mut_cov_matrix[, 2]
  mut_af[is.na(mut_af)] <- 0
  max_mut_af <- max(mut_af)

  # Set quantile cutoff equal to the observed mutant allele frequency
  q <- ceiling(mut_cov_matrix[, 2] * max_mut_af) - 1
  q[q < 0] <- 0
  if (all(q == 0)) {
    p <- 1
  } else {
    n <- sum(mut_cov_matrix[, 1])

    # Set multinomial probability from sequencing coverage
    prob <- mut_cov_matrix[, 2] / sum(mut_cov_matrix[, 2])
    prob[is.nan(prob)] <- 0
    p <- pmultinom(q, n, prob, lower.tail = FALSE)

    if (p < 0) {
      p <- 0
    } else if (p > 1) {
      p <- 1
    }
  }
  p
}

#' Merge multinucleotide indels (p < 0.05) into single indel events
#' Deletions are renamed to correct VCF style format, but insertions will need
#' to be renamed manually.
#' @param test_table data.frame; A table of variant information, rows are genome
#'                  positions
#' @return data.frame; test_table with indels merged
#' @export
merge_significant_indels <- function(test_table) {
  # Find an indel
  indel_idx <- find_indel(test_table)

  if (!is.null(indel_idx)) {
    # Recursively call merge_significant_indels, until all indels merged
    updated_table <- merge_indel(test_table, indel_idx)
    merge_significant_indels(updated_table)
  } else {
    test_table
  }
}

#' RSamtools pileup labels indels with "-". This function updates deletion.
#' fields in test_table with VCF formatting. This allows easy export to .vcf.
#' This function assumes all indels are small. Call merge_significant_indels
#' first to handle multinucleotide indels. Unfortunately Rsamtools doesn't
#' specify inserted bases so these have to be renamed manually.
#' @param test_table data.frame; A table of variant information, rows are genome
#'                  positions
#' @return data.frame
#' @export
rename_small_deletions <- function(test_table) {
  # Get indices of all deletions in table
  del_idx <- which(test_table$class == "deletion")
  ndel <- length(del_idx)
  # Output table
  new_t <- test_table

  for (i in 1:ndel) {
    j <- del_idx[i]
    if (test_table[j, "pos"] == 0) {
      prev_idx <- nrow[test_table]
    } else {
      prev_idx <- j - 1
      # If the first base is a small deletion, the new position is the very
      # last base
      if (prev_idx < 1) {
        prev_idx <- max(test_table$pos)
      }
    }
    new_t[j, "pos"] <- test_table[prev_idx, "pos"]
    new_t[j, "ref"] <- paste0(test_table[prev_idx, "ref"], test_table[j, "ref"])
    new_t[j, "alt"] <- test_table[j, "ref"]
  }

  new_t
}

# ==============================================================================
# END EXPORTS
# ==============================================================================

#' Generate a table of minor allele frequencies for each genome position in a
#' sample
#' @param pile matrix; allele counts at each genome position
#' @param consensus The consensus sequence for the samples
#' @return Numeric vector of allele frequencies, named by allele.
#' @importFrom fen.R.util ulapply
get_minor_allele_frequencies <- function(pile, consensus) {
  bp <- nrow(pile)
  minor_alleles <- get_minor_alleles(pile, consensus)
  mafs <- ulapply(1:bp, function(i) {
                    compute_allele_frequency(pile[i, ], minor_alleles[i])})
  mafs[is.nan(mafs)] <- 0
  mafs
}

#' Get the highest frequency allele in a vector of allele counts
#' WARNING: Returns the first value if tied
#' @param counts numeric vector; allele counts for a single genome position and
#'        sample
#' @return Character; One of {A, C, G, T, +, -}
get_highest_frequency_allele <- function(counts) {
  n <- length(counts)
  counts <- unlist(counts)
  if (n == 12) {
    counts <- destrand(counts)
  }

  # Get major allele
  allele <- names(counts)[which.max(counts)]
  allele
}

#' Get the allele with the second highest frequency for a single sample
#' WARNING: Returns the first value if tied
#' @param counts Numeric vector; Allele counts for a single genome position and
#'        sample
#' @param major_allele The major allele at this genome position
#' @return Character; One of {A, C, G, T, +, -}
get_minor_allele <- function(counts, major_allele) {
  # Minor allele is computed by finding the major allele, excluding the major
  # allele.
  non_consensus_alleles <- counts[! names(counts) == major_allele]
  get_highest_frequency_allele(non_consensus_alleles)
}

#' Calculate an allele's frequency from a vector of allele counts
#' @param counts Numeric vector; Allele counts for a single genome position and
#'        sample
#' @param allele Character; The target allele; One of {A, C, G, T, +, -}
#' @return Numeric
compute_allele_frequency <- function(counts, allele) {
  counts[allele] / sum(counts)
}

#' Calculate strand bias for a single genome position using fisher exact test
#' @param counts Numeric vector; Allele counts for a single genome position and
#'        sample
#' @param wild_type_allele Character; The wild type allele for this position
#' @param mutant_allele Character; The wild type allele for this position
#' @return Numeric; A p-value
#' @importFrom fen.R.util p_to_q
compute_strand_bias <- function(counts, wild_type_allele,
                                mutant_allele) {
  cnames <- names(counts)
  mut_idx <- which(grepl(paste0(mutant_allele, "_"), cnames, fixed = TRUE))
  wt_idx <- which(grepl(paste0(wild_type_allele, "_"), cnames, fixed = TRUE))
  fisher_matrix <- matrix(c(counts[mut_idx], counts[wt_idx]),
                          byrow = TRUE, ncol = 2)
  p <- fisher.test(fisher_matrix, workspace = 2e8)$p.value

  if (p > 1) {
    p <- 1
  } else if (p < 0) {
    p <- 0
  }

  phred <- p_to_q(p)
  phred
}

#' Carry out a 1-sided binomial test at each genome position for a sample
#' @param mut_cov_counts Matrix with two columns, first column is mutant
#'                       allele counts, second column is total sequencing
#'                       coverage. Rows are genome positions
#' @param seq_error_rate Numeric; Sequencing error rate estimate
#' @return Vector of p-values
binom_tests <- function(mut_cov_counts, seq_error_rate) {
  apply(mut_cov_counts, 1, function(mc) {
          binom.test(mc, p = seq_error_rate, alternative = "greater")$p.value
  })
}

#' Return TRUE if only one of a list of p-values is < 0.05
is_unique <- function(ps) {
  sorted <- sort(ps)
  (sorted[1] < 0.05) && (sorted[2] > 0.05)
}

#' Split a vector by its runs of consecutive values
#' E.g partition_by_runs(c(1, 2, 4, 5)) == list(c(1, 2), c(4, 5))
#' @param x Numeric vector
#' @return List of numeric vectors
partition_by_runs <- function(x) {
  split(x, cumsum(seq_along(x) %in% (which(diff(x) > 1) + 1)))
}

#' Find runs of consecutive values with length > 1
#' @param x Numeric vector
#' @return List of numeric vectors
find_runs <- function(x) {
  # Partition idx by runs
  par <- partition_by_runs(x)
  # Return runs > 1
  runs <- par[lengths(par) > 1]
  runs
}

#' Find a multinucleotide indel in a table of variants
#' @param test_table data.frame; A table of variant information, rows are genome
#'                  positions
#' @return Numeric vector; row indices of a multinucleotide indel. NULL if no
#'         indels can be found.
find_indel <- function(test_table) {

  # First look for insertions
  is_ins <- test_table[, "alt"] == "+"
  is_sig <- test_table$p_value < 0.05

  insertion_idx <- which(is_ins & is_sig)
  insertion_run <- find_runs(insertion_idx)

  if (length(insertion_run) == 0) {
    # If we've found all insertions, look for deletions
    is_del <- test_table[, "alt"] == "-"
    deletion_idx <- which(is_del & is_sig)
    deletion_run <- find_runs(deletion_idx)

    # If all indels have been merged, return NULL
    if (length(deletion_run) == 0) {
      indel_idx <- NULL
    } else {
      indel_idx <- deletion_run[[1]]
    }

  } else {
    indel_idx <- insertion_run[[1]]
  }

  indel_idx
}

#' Merge rows of test_table given by indel_idx into a single row.
#' @param test_table data.frame; A table of variant information, rows are genome
#'                  positions
#' @param indel_idx Numeric vector (length > 1)
#' @return A test_table
merge_indel <- function(test_table, indel_idx) {
  pre_row <- test_table[indel_idx[1] - 1, ]
  indel_row <- test_table[indel_idx,]

  # Determine alternative and reference alleles based upon indel class
  in_or_del <- unique(indel_row$alt)
  pre_ref <- pre_row$ref
  if (in_or_del == "-") {
    class <- "deletion"
    del_ref <- paste0(indel_row$ref, collapse = "")
    ref <- paste0(pre_ref, del_ref, collapse = "")
    alt <- pre_ref
  } else if (in_or_del == "+") {
    class <- "insertion"
    ref <- pre_ref
    alt <- paste0(pre_ref, "+", collapse = "")
  } else {
    stop("These indices don't correspond to an indel.")
  }

  # The new merged row.
  new_row <- list(
                  # species
                  unique(as.character(indel_row$species)),
                  # population
                  unique(as.character(indel_row$population)),
                  # genotype
                  unique(as.character(indel_row$genotype)),
                  # mutant sample
                  unique(as.character(indel_row$sample))[1],
                  # position
                  pre_row$pos,
                  # coverage
                  mean(indel_row$coverage),
                  # reference allele
                  ref,
                  # variant allele
                  alt,
                  # mutation class
                  class,
                  # mutant allele frequency
                  mean(indel_row$af),
                  # mutant sample allele frequency - mean allele frequency
                  mean(indel_row$af_diff),
                  # strand bias
                  mean(indel_row$strand_bias),
                  # is variant unique?
                  all(indel_row$unique),
                  # is variant low coverage?
                  mean(indel_row$coverage_proportion),
                  # p value
                  mean(indel_row$p_value))
  new_test_table <- replace_rows(test_table, indel_idx, new_row)
  new_test_table
}
