#!/usr/bin/env Rscript
"Produce a spliced pileup from two .bam files resulting from alignment to
original and rotated mtDNA reference sequences.

Usage:
  produce_spliced_pileups.R --output PILE OGBAM ROTBAM

Options:
  -o --output     Output file" -> doc

library(docopt)

main <- function(og, rot, out) {
  pile <- megadaph.mtdna::construct_spliced_pileup(og, rot,
                                                   distinguish_strands = TRUE)
  pile
}

if (!interactive()) {
  opts <- docopt::docopt(doc)
  library(megadaph.mtdna)
  og <- unlist(opts["OGBAM"])
  rot <- unlist(opts["ROTBAM"])
  out <- unlist(opts["PILE"])
  pile <- main(og, rot, out)
  write.csv(pile, out, quote = FALSE, row.names = FALSE)
}
