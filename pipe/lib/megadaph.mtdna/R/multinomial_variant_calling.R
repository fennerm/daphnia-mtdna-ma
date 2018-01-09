# ==============================================================================
# EXPORTS
# ==============================================================================

#' Splice the original and rotated BAM files into a single allele count table
#' @param og_bam .bam files from alignment to the unrotated reference
#'                sequences
#' @param rot_bam .bam files from alignment to the rotated reference
#'                sequences
#' @param min_base_quality Numeric; minimum Phred base quality for base to be
#'        included in table
#' @param distinguish_strands; If TRUE base counts on the +/- strands are
#'        counted separately
#' @return A data.frame
#' @export
construct_spliced_pileup <- function(og_bam, rot_bam, min_base_quality=30,
                                     distinguish_strands=FALSE) {
  # Create individual pileups
  og_pile <- construct_pileup(og_bam, min_base_quality, distinguish_strands)
  rot_pile <- construct_pileup(rot_bam, min_base_quality, distinguish_strands)

  # Splice the pileups
  spliced_pile <- splice_pileups(og_pile, rot_pile, species)
  spliced_pile
}

#' Create an allele count pileup with RSamtools
#' @param bam Path to a .bam file
#' @param min_base_quality All read positions with phred score <
#'        min_base_quality will be excluded
#' @param distinguish_strands; If TRUE, base counts on the +/- strands are
#'        counted separately
#' @param min_nucleotide_depth; Minimum allele count to be included in pileup
#' @return A data.table
#' @importFrom Rsamtools PileupParam pileup
#' @importFrom data.table setDT
#' @importFrom reshape2 dcast
#' @export
construct_pileup <- function(bam, min_base_quality = 30,
                             distinguish_strands = FALSE,
                             min_nucleotide_depth = 1) {
  pileup_param <- PileupParam(max_depth = 1000000,
                              distinguish_strands = distinguish_strands,
                              distinguish_nucleotides = TRUE,
                              ignore_query_Ns = TRUE,
                              min_nucleotide_depth = min_nucleotide_depth,
                              include_deletions = TRUE,
                              include_insertions = TRUE,
                              min_base_quality = min_base_quality)

  pile <- pileup(bam, pileupParam = pileup_param)

  # Cast the table to wide format
  if (distinguish_strands) {
    pile_wide <- dcast(pile,
                       seqnames + pos ~ nucleotide + strand,
                       value.var = "count",
                       fill = 0L)
  } else {
    pile_wide <- dcast(pile,
                       seqnames + pos ~ nucleotide,
                       value.var = "count",
                       fill = 0L)
  }
  pile_wide <- setDT(pile_wide)
  # Remove unnecessary columns
  pile_wide[, c("seqnames", "pos") := NULL]
  pile_wide
}

#' Create a consensus sequence from a set of nucleotide count pileups
#' @param piles List of data.table; Each list item contains allele counts at
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
#' @param piles List of data.table; Each list item contains allele counts at
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
                                v <- ulapply(mafs, "[", i)
                                vnames <- lapply(mafs, names)
                                vnames[[which.max(v)]][i]
                       })
  mutant_consensus
}

#' Calculate mutant allele frequency at each genome position for a sample
#' @param pile A data.table; allele counts at each genome position
#' @param mutant_consensus; An alternate allele consensus sequence produced by
#'        create_mutant_consensus
#' @return Numeric vector
#' @importFrom fen.R.util ulapply
#' @export
compute_mutant_allele_frequencies <- function(pile, mutant_consensus) {
  bp <- nrow(pile)
  ulapply(1:bp, function(i) {
            counts <- unlist(pile[i, ])
            compute_allele_frequency(counts, mutant_consensus[i])
                       })
}

#' Calculate strand bias for each genome position using fisher tests
#' @param stranded_piles data.table; Allele count pileups separated by strand
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
#' @param pile A data.table; allele counts at each genome position
#' @param mutant_consensus; An alternate allele consensus sequence produced by
#'        create_mutant_consensus
#' @return An nx2 matrix, where columns are mutant and wild-type allele counts,
#'         and rows are genome positions.
#' @export
convert_to_mut_wt_counts <- function(pile, mutant_consensus) {
  mut_counts <- get_allele_counts(pile, mutant_consensus)
  cov <- get_coverage(pile)
  wt_count <- cov - mut_counts
  cbind(mut_counts, non_mut_counts)
}

#' Get the allele counts for a set alleles
#' @param pile A data.table; allele counts at each genome position
#' @param alleles Character vector
#' @return Numeric vector
#' @importFrom fen.R.util ulapply
#' @export
get_allele_counts <- function(pile, alleles) {
  bp <- nrow(pile)
  counts <- ulapply(1:bp, function(i) {
                      pile[i, alleles[i], with = FALSE]
                                    })
  counts

}

#' Get the major allele at each genome position for a single sample
#' WARNING: Returns the first value if tied
#' @param pile data.table; allele counts at each genome position
#' @return Character; Vector of major alleles, each allele one of
#'         {A, C, G, T, +, -}
#' @export
get_major_alleles <- function(pile) {
  apply(pile, 1, get_highest_frequency_allele)
}

#' Get the alleles with the second highest frequency at each position
#' @param pile data.table; allele counts at each genome position
#' @param consensus The consensus sequence for the samples
#' @return Character; Vector of minor alleles, each allele one of
#'         {A, C, G, T, +, -}
#' @importFrom fen.R.util ulapply
#' @export
get_minor_alleles <- function(pile, consensus) {
  bp <- nrow(pile)
  ulapply(1:bp, function(i) {
            get_minor_allele(unlist(pile[i, ]),
                             consensus[i])
  })
}

#' Get sequencing coverage for each genome position
#' @param pile A data.table; allele counts at each genome position
#' @return A numeric vector
#' @export
get_coverage <- function(pile) {
  # Function for getting coverage of a single sample.
  get_sample_coverage <- function(p) {
    unlist(apply(pile, 1, sum))
  }

  # If list of data.tables, apply to each.
  if (is.null(ncol(pile))) {
    lapply(pile, get_sample_coverage)
  } else {
    get_sample_coverage(pile)
  }
}

#' Convert an allele count pileup to a table of mutant vs. coverage
#' @param pile A data.table; allele counts at each genome position
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
#' @param mut_cov_counts List of data.frames; Each data.frame has two columns,
#'                       first column is mutant allele counts, second column is
#'                       total sequencing coverage. Rows are genome positions
#' @param seq_error_rates Numeric vector; Sequencing error rate estimates for
#'                        each sample
#' @return Logical vector
determine_unique <- function(mut_cov_counts, seq_error_rates) {
  nsamples <- length(mut_cov_counts)

  # Run binomial test on each genome position for each sample
  p <- mapply(binom_tests, mut_cov_counts, seq_error_rates, SIMPLIFY = FALSE)

  # Convert list to multisample matrix
  p_matrix <- matrix(unlist(p), ncol=nsamples)

  unique <- apply(p_matrix, 1, is_unique)
  unique
}

# ==============================================================================
# END EXPORTS
# ==============================================================================

#' Splice the 'original' and 'rotated' pileups into a single pileup
#' @param og_pile Pileup from alignment to the 'original' reference sequence
#' @param rot_pile Pileup from alignment to the 'rotated' reference sequence
#' @return A spliced pileup containing the middle halves of both the 'rotated'
#' and 'original' pileups
#' @importFrom data.table rbindlist
splice_pileups <- function(og_pile, rot_pile, rot_ref) {
  bp <- nrow(og_pile)
  splx <- compute_split_indices(bp)

  data_type <- class(og_pile)[1]

  if (data_type == "matrix") {
    spliced <- rbind(rot_pile[splx[[1]], ], og_pile[splx[[2]],],
                     rot_pile[splx[[3]],])
  } else if (data_type == "data.table") {
    spliced <- rbindlist(list(rot_pile[splx[[1]],], og_pile[splx[[2]],],
                              rot_pile[splx[[3]],]))
  } else if (data_type %in% c("numeric", "integer")) {
    spliced <- c(rot_pile[splx[[1]]], og_pile[splx[[2]]], rot_pile[splx[[3]]])
  }
  spliced
}

#' Compute the indices for splicing the 'original' and 'rotated' pileups
#' @param sequence_length Integer; The length of the reference sequence
#' @return List, length=3; First and third list items are indices from 'rotated'
#' pileup. Second item is indices into 'original' pileup.
compute_split_indices <- function(sequence_length) {
  spl <- floor(seq(1, sequence_length, len = 5))
  rot1 <- og_to_rot(spl[1], sequence_length):og_to_rot(spl[2] - 1,
                                                       sequence_length)
  og <- spl[2]:(spl[4] - 1)
  rot2 <- og_to_rot(spl[4], sequence_length):og_to_rot(spl[5], sequence_length)
  list(rot1, og, rot2)
}

#' Generate a table of minor allele frequencies for each genome position in a
#' sample
#' @param pile data.table; allele counts at each genome position
#' @param consensus The consensus sequence for the samples
#' @return Numeric vector of allele frequencies, named by allele.
#' @importFrom fen.R.util ulapply
get_minor_allele_frequencies <- function(pile, consensus) {
  bp <- nrow(pile)
  minor_alleles <- get_minor_alleles(pile, consensus)
  mafs <- ulapply(1:bp, function(i) {
                    compute_allele_frequency(unlist(pile[i, ]),
                                             minor_alleles[i])})
  mafs
}

#' Get the highest frequency allele in a vector of allele counts
#' WARNING: Returns the first value if tied
#' @param counts Numeric vector; Allele counts for a single genome position and
#'        sample
#' @return Character; One of {A, C, G, T, +, -)
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

#' Consolidate strand counts in a strand-split allele count table
#' @param counts data.table; Table with nucleotide and indel counts for each
#'        genome position split by strand
#' @return data.table; A allele count table with the stranded counts summed
#' @export
#' @importFrom data.table setDT
destrand <- function(counts) {
  idx1 <- seq(1, 11, by = 2)
  idx2 <- idx1 + 1

  if (is.numeric(counts)) {
    destranded <- counts[idx1] + counts[idx2]
  } else {
    sum_counts <- function(x1, x2) {
      counts[, x1, with = FALSE] + counts[, x2, with = FALSE]
    }
    destranded <- mapply(sum_counts, idx1, idx2)
    destranded <- setDT(destranded)
  }
  names(destranded) <- c("A", "C", "G", "T", "-", "+")
  destranded
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
#' @param mut_cov_counts Data.frame with two columns, first column is mutant
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


# Merge rows of test_table given by indel_idx into a single row.
# Given:
#   test_table - A table generated by construct_test_tables
#   indel_idx - A numeric vector (length > 1)
# Return:
#   A test_table with indel_idx rows merged into single row.
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
    alt <- paste0(pre_ref, "+", collapse = "")
  } else {
    stop("These indices don't correspond to an indel.")
  }

  # The new merged row.
  new_row <- list(
                  # species
                  unique(as.character(indel_row$species)),
                  # genotype
                  unique(as.character(indel_row$genotype)),
                  # isolate
                  unique(as.character(indel_row$isolate)),
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


#' @export
#' @importFrom qvalue qvalue
mult_comparisons_correct <- function(test_tables, fdr_level=0.05) {
  merged_table <<- do.call(rbind, test_tables)
  # filtered_table <- merged_table[which(!merged_table$low_coverage),]
  fdr <- qvalue(merged_table$p_value, fdr.level=fdr_level)
  significant <- fdr$significant
  q_value <- fdr$qvalues
  fdr_table <- cbind(merged_table, q_value, significant)
  fdr_table <- data.frame(fdr_table, stringsAsFactors=FALSE)
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

calc_coverages <- function(test_tables) {
  coverages <- sapply(test_tables, function(x) mean(x$coverage))
  coverages
}

partition_by_runs <- function(x) {
  split(x, cumsum(seq_along(x) %in%
                  (which(diff(x) > 1) + 1)))
}

# Find runs of identical values
# Given:
#   idx - A numeric vector
# Return:
#   A list of numeric vectors of values in idx. Each list item is a run.
find_runs <- function(x) {
  # Partition idx by runs
  par <- partition_by_runs(x)
  # Return runs > 1
  runs <- par[lengths(par) > 1]
  runs
}

# Find a multinucleotide indel in a test table
# Given:
#   test_table - A table generated by construct_test_tables
# Return:
#   The indices of an arbitrary multinucleotide indel in test_table. If no
#   indels found, NULL is returned.
#' @export
find_indel <- function(test_table) {

  # First look for insertions
  is_ins <- test_table[, "alt"] == "+"
  is_sig <- test_table$p_value < 0.05

  insertion_idx <- which(is_ins & is_sig)
  insertion_run <- find_runs(insertion_idx)

  if (length(insertion_run) == 0) {
    is_del <- test_table[, "alt"] == "+"

    # If we've found all insertions, look for deletions
    deletion_idx <- which(is_del & is_sig)
    deletion_run <- find_runs(deletion_idx)

    # If all indels have been merged, return NULL
    if (length(deletion_run) == 0) {
      NULL
    } else {
      deletion_run[[1]]
    }
  } else {
    insertion_run[[1]]
  }
}

## Find runs of significant (p < 0.05) insertions or deletions and merge them
## into single indel events. Deletions are renamed to correct VCF style format,
## but insertions will need to be renamed manually.
## Given:
##   test_table - A table generated by construct_test_tables
## Return:
##   A test_table
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

## RSamtools pileup labels indels with "-". This function updates deletion.
## fields in test_table with VCF formatting. This allows easy export to .vcf.
## This function assumes all indels are small. Call merge_significant_indels
## first to handle multinucleotide indels. Unfortunately Rsamtools doesn't
## specify inserted bases so these have to be renamed manually.
## Given:
##   test_table - A table generated by construct_test_tables
## Return:
##   A test_table
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
    }
    new_t[j, "pos"] <- test_table[prev_idx, "pos"]
    new_t[j, "ref"] <- paste0(test_table[prev_idx, "ref"],
                              test_table[j, "ref"])
    new_t[j, "alt"] <- test_table[j, "ref"]
  }

  new_t
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


#' Convert allele count pileups to mutant vs wildtype pileups
#'
#' @param piles list of data.table; each list item contains allele counts at
#'        each genome position for a single sample.
#' @param mutant_consensus; An alternate allele consensus sequence produced by
#'        `create_mutant_consensus`
#' @export
piles_to_mut_wt <- function(piles, mutant_consensus) {
  bp <- length(mutant_consensus)
  mut_wt_counts <- lapply(piles, convert_to_mut_wt_counts, mutant_consensus)
  mut_wt_matrices <- lapply(1:bp, function(i) {
                              mat <- sapply(mut_wt_counts, function(x) {
                                              x[i, ]
                  })
                              t(mat)
                              })
  mut_wt_matrices
}

#' Given a column index into an unstranded pileup, return the equiv. stranded
#' indices
.stranded_indices <- function(idx) {
  str_idx <- list(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 10), c(11, 12))
  str_idx <- matrix(unlist(str_idx[idx]), ncol = 2, byrow = TRUE)
  str_idx
}
