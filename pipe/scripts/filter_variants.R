#!/usr/bin/env Rscript
"Exclude variants which do not meet inclusion criteria from the test tables

The test tables are expected to include q-values from multiple comparisons
correction.

Usage:
  filter_variants.R --max_strand_bias=FLOAT --output=CSV TEST_TABLE

Options:
  -s --max_strand_bias=FLOAT  Maximum strand bias allowed for included variants.
  -o --output=CSV    Output file name
" -> doc


main <- function(input_file, output_file, max_strand_bias) {
  fdr_table <- read.csv(input_file, stringsAsFactors = FALSE)
  var_table <- megadaph.mtdna::filter_variants(fdr_table, max_strand_bias)
  invisible(write.csv(var_table, output_file, row.names = FALSE, quote = FALSE))
}


if (!interactive()) {
  opts <- docopt::docopt(doc)
  max_strand_bias <- as.numeric(unlist(opts["max_strand_bias"]))
  input_file <- normalizePath(unlist(opts["TEST_TABLE"]))
  output_file <- unlist(opts["output"])
  main(input_file, output_file, max_strand_bias)
}
