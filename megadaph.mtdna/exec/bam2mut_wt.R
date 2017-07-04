#!/usr/bin/env Rscript
# pile_to_mut_wt <- function() {
# }
#
# bams_to_piles <- function() {
#
# }
#
# create piledir and mut_wt_dir
#
# read file to pile (all threads 1 pile)
# save piles
#
# create mut_consensus (all threads 1 isolate)
# at each position:
#     create mut_wt




# Parameters for creating pileup from bam files.
pileup_param <- Rsamtools::PileupParam(max_depth = 1000000,
                                       distinguish_strands = TRUE,
                                       distinguish_nucleotides = TRUE,
                                       ignore_query_Ns = TRUE,
                                       include_deletions = TRUE,
                                       include_insertions = TRUE)
