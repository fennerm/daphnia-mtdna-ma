#!/usr/bin/env Rscript
#' Variants were called using a novel contingency matrix approach.
#' For each group of samples with a common ancestor genotype:
#' 1. Determine the major and minor allele at every locus in each sample.
#'    - The major allele is defined as the nucleotide which is present in the
#'      the most reads at a site.
#' 2. Determine the consensus major allele at each site across all samples.
#' 3. Define the minor allele as the allele which reached the highest allele
#'    frequency in a single sample, other than the major allele.
#' 4. For each site in the genome, produce a contingency table of the form
#'              | Major | Minor |
#'     Sample 1 |   -   |   -   |
#'        ...   |   -   |   -   |
#'     Sample N |   -   |   -   |
#'
#'    where 'Major' and 'Minor' are the major and minor allele frequency in
#'    each sample.
#' 5. Test the null hypothesis that the allele frequencies at a particular site
#'    are equal across samples.
#'    - For this purpose we adapted a statistical test based on the expected
#'      distribution of the maximum of a hypergeometric distribution:
#'      The exact distribution of the maximum, minimum and the range of
#'      Multinomial/Dirichlet and Multivariate Hypergeometric frequencies,
#'      Corrado, 2011
library(megadaph.mtdna)
library(fen.R.util)

DELETE <- function() {
  library(fen.R.util)
  og_bams <- c("../../../../../data/alignments/magna/F/og/FA12.sorted.grouped.viterbi.sorted.realign.mrkdup.bam", "../../../../../data/alignments/magna/F/og/FA1.sorted.grouped.viterbi.sorted.realign.mrkdup.bam")
  rot_bams <- c("../../../../../data/alignments/magna/F/rot/FA12.sorted.grouped.viterbi.sorted.realign.mrkdup.bam", "../../../../../data/alignments/magna/F/rot/FA1.sorted.grouped.viterbi.sorted.realign.mrkdup.bam")
  og_bam <- og_bams[2]
  rot_bam <- rot_bams[2]
  piles <- mapply(construct_spliced_pileup, og_bams, rot_bams, SIMPLIFY = FALSE)
  pile <- piles[[1]]
}

main <- function(og_bams, rot_bams, rot_ref, seq_err) {
  nsamples <- length(og_bams)

  cat("Building allele count pileups \n")
  piles <- mapply(construct_spliced_pileup, og_bams, rot_bams, SIMPLIFY = FALSE)

  cat("Building consensus sequence \n")
  consensus <- create_consensus(piles)
  bp <- length(consensus)

  consensus_string <- paste0(consensus, collapse = "")
  fileconn <- file(fasta)
  writeLines(consensus_string, fileconn)
  close(fileconn)

  cat("Determining possible mutant alleles at each position \n")
  mutant_consensus <- create_mutant_consensus(piles, consensus)

  cat("Calculating allele count statistics \n")
  wt_counts <- lapply(piles, get_allele_counts, consensus)
  mutant_counts <- lapply(piles, get_allele_counts, mutant_consensus)
  coverage <- lapply(piles, get_coverage)
  mean_coverage <- unlist(apply_across(coverage, mean))
  mut_cov_counts <- mapply(cbind, mutant_counts, coverage, SIMPLIFY = FALSE)
  overall_coverage <- mean(mean_coverage)
  coverage_proportion <- mean_coverage / overall_coverage

  cat("Calculating mutant allele frequency for each site in each sample \n")
  mut_afs <- mapply("/", mutant_counts, coverage)

  cat("Determining which sample has the highest mutant allele frequency at each
      site")
  mutant_samples <- apply(mut_afs, 1, which.max)

  cat("Calculating strand bias \n")
  stranded_piles <- mapply(construct_spliced_pileup, og_bams, rot_bams,
                           distinguish_strands = TRUE, SIMPLIFY = FALSE)
  strand_bias <- compute_all_strand_bias(stranded_piles, consensus,
                                         mutant_consensus, mutant_samples)

  cat("Determining unique vs. shared variants \n")
  unique <- determine_unique(mut_cov_counts, seq_err)

  cat("Gathering sample metadata \n")
  samples <- sapply(og_bams, get_sample)
  isolate <- unique(get_isolate(og_bams))
  genotype <- unique(sapply(og_bams, get_genotype))
  species <- unique(get_species(og_bams))

  ## HERE
  var_samples <- ulapply(c(1:bp), function(i) samples[max_sample[i]])

  # Get variant allele frequencies
  var_afs <- sapply(1:bp, function(i) {
                      mut_afs[i, max_sample[i]]
                                         })

  # Get difference between mutant sample and non-mutant sample allele
  # frequencies
  diff_afs <- sapply(1:bp, function(i) {
                       var_afs[i] - mean(mut_afs[i,-max_sample[i]])
                                         })

  # Get class of mutation
  mut_class <- sapply(mutant_consensus, function(allele) {
                        if (allele == "-") {
                          "deletion"
                        } else if (allele == "+") {
                          "insertion"
                        } else {
                          "snv"
                        }
  mut_wt_matrices <- lapply(1:bp, function(i) {
                              mat <- sapply(mut_wt_counts, function(x) {
                                              x[i,]
                            })
                              t(mat)
                                         })
                                             })

  # Return output table
  test_table <- cbind(as.data.frame(rep(species, bp), stringsAsFactors=FALSE),
                      rep(genotype, bp), isolate, var_samples, c(1:bp),
                      mean_coverage, consensus, mutant_consensus,
                      mut_class, var_afs, diff_afs,
                      sb, unique, coverage_proportion,
                      stringsAsFactors=FALSE)
  colnames(test_table) <- c("species", "genotype", "isolate", "sample", "pos",
                            "coverage", "ref", "alt", "class", "af", "af_diff",
                            "strand_bias", "unique", "coverage_proportion")
  list(piles, mut_wt_matrices, test_table)


  tabs <<- lapply(1:nsamples, function(i) {
                    construct_tables(og_bams[[i]], rot_bams[[i]],
                                     seq_err_rates[[i]], subsample)}
  )
  piles <<- lapply(tabs, "[[", 1)
  mut_wt_matrices <<- lapply(tabs, "[[", 2)
  test_tables <<- lapply(tabs, "[[", 3)
  merged_table <- do.call(rbind, test_tables)
  p_value <- readRDS("../tables/corrado_p.Rds")
  merged_table <- cbind(merged_table, p_value)
  test_tables <- split(merged_table, merged_table$isolate)
  formatted_test_tables <- lapply(test_tables, merge_significant_indels)
  formatted_test_tables <<- lapply(formatted_test_tables, rename_small_deletions)
  fdr_table <<- mult_comparisons_correct(formatted_test_tables,
                                         fdr_level=0.01)
  var_table <<- get_stat_significant(fdr_table)
  var_table
}

if (interactive()) {

}
# args = commandArgs(trailingOnly=TRUE)
# mut_wt_matrices <- readRDS(as.character(args[1]))
# threads <- as.numeric(args[2])
# sq <- 1:length(mut_wt_matrices)
# partitions <- split(sq, cut(sq, 40, labels=FALSE))
# cl <- makeCluster(threads, outfile="cluster.log")
# clusterExport(cl, c("logmatrix_mult", "independent_assumption_test",
#                     "normalize", "cull", "build_stochastic_matrix",
#                     "build_choose_matrix", "probability_proportion",
#                     "log_sum_exp", "corrado_test.R", "mut_wt_matrices",
#                     "partitions"))
# p <- unlist(parLapply(cl, partitions, function(is) {
#     unlist(lapply(mut_wt_matrices[is], corrado_test))
# }))
# saveRDS(p, "corrado_p.Rds")
# stopCluster(cl)
