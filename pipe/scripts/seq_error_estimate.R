#!/usr/bin/env Rscript
"Calculate sequencing error rates for heterozygous samples.

Error rates are calculated from nuclear genome alignments so that heteroplasmy
does not have to be accounted for.

Usage: seq_error_estimate.R INPUT_BAM OUTPUT_CSV" -> doc

library(data.table)
library(docopt)
library(Rsamtools)

## Calculate the sequencing error rate for a single bam file.
## param: bam A path to a bam file
## return: Numeric
main <- function(bam) {
  pileup_param <- Rsamtools::PileupParam(max_depth = 1000000,
                                         distinguish_strands = FALSE,
                                         distinguish_nucleotides = TRUE,
                                         ignore_query_Ns = TRUE,
                                         include_deletions = TRUE,
                                         include_insertions = TRUE)
  # Create base count pilep and convert to a data.table
  pile <- data.table::setDT(Rsamtools::pileup(bam, pileupParam = pileup_param))

  # Convert from long to wide format and exclude position and c'some fields
  pile <- data.table::dcast(pile, seqnames + pos ~ nucleotide,
                            value.var="count", fill = 0L)
  pile <- pile[, 3:length(pile)]

  consensus <- get_major_alleles(pile)
  minor_alleles <- get_minor_alleles(pile, consensus)
  minor_allele_counts <- get_allele_counts(pile, minor_alleles)

  # Coverage at each position
  sum_counts <- apply(pile, 1, sum)

  # Homozygous positions
  # 0.2 chosen as cutoff since pbinom is approx 0.001 for depth of 20.
  cutoff <- 0.2 * sum_counts
  homo <- which((minor_allele_counts < cutoff) & (sum_counts > 39))

  # Sequencing error estimate
  seq_err <- mean(minor_allele_counts[homo] / sum_counts[homo])
  seq_err
}

if (!interactive()) {
  opts <- docopt::docopt(doc)
  sample <- megadaph.mtdna::get_sample(opts["INPUT_BAM"])
  seq_err_rate <- main(opts["INPUT_BAM"])
  seq_err_table <- data.frame(sample = sample, sequencing_error = seq_err_rate,
                              stringsAsFactors = FALSE)
  write.csv(seq_err_rate, opts["OUTPUT_CSV"], quote=FALSE, row.names = FALSE)
}
