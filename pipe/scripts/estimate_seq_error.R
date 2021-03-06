#!/usr/bin/env Rscript
"Calculate sequencing error rates for heterozygous samples.

Error rates are calculated from nuclear genome alignments so that heteroplasmy
does not have to be accounted for.

Usage: seq_error_estimate.R INPUT_BAM OUTPUT_CSV" -> doc

main <- function(bam) {
  pileup_param <- Rsamtools::PileupParam(max_depth = 1000000,
    min_base_quality = 30,
    distinguish_strands = FALSE,
    distinguish_nucleotides = TRUE,
    ignore_query_Ns = TRUE,
    include_deletions = TRUE,
    include_insertions = TRUE)
  cat("Creating pileup \n")
  pile <- Rsamtools::pileup(bam, pileupParam = pileup_param)
  if (nrow(pile) == 0) {
    seq_err <- 0
    stdev <- 0
  } else {
    pile <- data.table::setDT(pile)

    cat("Converting from long to wide format \n")
    pile <- data.table::dcast(pile, seqnames + pos ~ nucleotide,
      value.var="count", fill = 0L)

    cat("Excluding extra columns \n")
    pile <- pile[, 3:length(pile)]
    pile <- as.matrix(pile)

    cat("Determining major alleles \n")
    consensus <- megadaph.mtdna::get_major_alleles(pile)

    cat("Determining minor alleles \n")
    minor_alleles <- megadaph.mtdna::get_minor_alleles(pile, consensus)
    minor_allele_counts <- megadaph.mtdna::get_allele_counts(pile, minor_alleles)

    # Coverage at each position
    sum_counts <- apply(pile, 1, sum)

    cat("Determining homozygous positions \n")
    # Homozygous positions
    # 0.2 chosen as cutoff since pbinom is approx 0.001 for depth of 20.
    cutoff <- 0.2 * sum_counts

    homo <- which((minor_allele_counts < cutoff) & (sum_counts > 39) &
      (sum_counts < 200))

    if (length(homo) > 0) {
      # Sequencing error estimate
      minor_fraction <- minor_allele_counts[homo] / sum_counts[homo]
      cat("Calculating sequencing error \n")
      seq_err <- mean(minor_fraction)
      stdev <- sd(minor_fraction)
    } else {
      # If no homozygous positions found, set sequencing error to 0
      # This is primarily to ensure that a value is still returned for small test
      # files.
      seq_err <- 0
      stdev <- 0
    }
  }

  c(seq_err, stdev)
}

if (!interactive()) {
  opts <- docopt::docopt(doc)
  input_bam <- unlist(opts["INPUT_BAM"])
  output_csv <- unlist(opts["OUTPUT_CSV"])
  smp <- megadaph.mtdna::get_sample(input_bam)
  seq_err_stats <- main(input_bam)

  cat("Writing output \n")
  seq_err_table <- data.frame(sample = smp,
    sequencing_error = seq_err_stats[1],
    stdev = seq_err_stats[2],
    stringsAsFactors = FALSE)
  write.csv(seq_err_table, output_csv, quote = FALSE, row.names = FALSE)
}
