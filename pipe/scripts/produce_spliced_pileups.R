#!/usr/bin/env Rscript
"Produce a spliced pileup from two .bam files resulting from alignment to
original and rotated mtDNA reference sequences.

Usage:
  produce_spliced_pileups.R --output=PILE --bp=N OGBAM ROTBAM

Options:
  -o --output=PILE     Output file
  -n --bp=N            Number of bases in reference sequence
" -> doc

library(docopt)

library(megadaph.mtdna)

main <- function(og, rot, bp) {
  pile <- megadaph.mtdna::construct_spliced_pileup(og, rot, bp,
                                                   distinguish_strands = TRUE)
  pile
}

if (!interactive()) {
  opts <- docopt::docopt(doc)
  bp <- as.numeric(opts["bp"])
  pile <- main(opts["OGBAM"], opts["ROTBAM"], bp)
  write.csv(pile, out, quote = FALSE, row.names = FALSE)
}
