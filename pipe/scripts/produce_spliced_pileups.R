#!/usr/bin/env Rscript
"Produce a spliced pileup from two .bam files resulting from alignment to
original and rotated mtDNA reference sequences.

Usage:
  produce_spliced_pileups.R --output PILE --bp N OGBAM ROTBAM

Options:
  -o --output     Output file
  -n --bp         Number of bases in reference sequence
" -> doc

library(docopt)
opts <- docopt::docopt(doc)

library(megadaph.mtdna)

main <- function(og, rot, bp, out) {
  pile <- megadaph.mtdna::construct_spliced_pileup(og, rot, bp,
                                                   distinguish_strands = TRUE)
  pile
}

if (!interactive()) {
  og <- unlist(opts["OGBAM"])
  rot <- unlist(opts["ROTBAM"])
  out <- unlist(opts["--output"])
  bp <- int(unlist(opts["--bp"]))
  pile <- main(og, rot, out)
  write.csv(pile, out, quote = FALSE, row.names = FALSE)
}
