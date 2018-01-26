#!/usr/bin/env Rscript
"Add q values to a test table

Usage:
  correct_for_multiple_comparisons.R [--fdr=FDR_LEVEL] --outdir=DIR TEST_TABLE ...

Options:
  -o --outdir=DIR    Output directory
  -d --fdr=FDR_LEVEL False discovery rate level [default: 0.05]
" -> doc


write_output <- function(table, output_name, outdir) {
  write.csv(table, file = file.path(outdir, output_name), row.names = FALSE,
            quote = FALSE)
}

main <- function(input_files, fdr_level, outdir) {
  output_names <- lapply(input_files, basename)
  # Read input files
  test_tables <- lapply(input_files, read.csv, stringsAsFactors = FALSE)
  # Add q-values to the tables
  fdr_tables <- megadaph.mtdna::mult_comparisons_correct(test_tables, fdr_level)
  # Write output files
  invisible(mapply(write_output, fdr_tables, output_names, outdir = outdir))
}


if (!interactive()) {
  opts <- docopt::docopt(doc)
  fdr <- as.numeric(opts["fdr"])
  test_tables <- sapply(opts["TEST_TABLE"], normalizePath)
  outdir <- normalizePath(unlist(opts["outdir"]))
  main(test_tables, fdr, outdir)
}
