#!/usr/bin/env Rscript
"Merge multiple single sample sequencing error files from a single genotype

Usage: merge_sequencing_error_files.R --output OUTPUT CSV ...

Options:
  -o --output     Output file" -> doc

main <- function(seq_err_files) {
  seq_err <- lapply(seq_err_files, read.csv, stringsAsFactors = FALSE)
  merged_table <- do.call(rbind, seq_err)
  merged_table <- merged_table[order(merged_table$sample),]
  merged_table
}

if (!interactive()) {
  opts <- docopt::docopt(doc)
  seq_err_files <- unlist(opts["CSV"])
  output <- unlist(opts["OUTPUT"])
  merged_table <- main(seq_err_files)
  write.csv(merged_table, file = output, row.names = FALSE, quote = FALSE)
}
