#!/usr/bin/env Rscript
"Annotate variants with predicted effects and ts/tv

Usage:
annotate_variants.R --output FILE --config=CONF VAR_TABLE

Options:
  -c --config=CONF  snpEff config file
  -o --output=FILE  Output csv file
" -> doc
library(tidyverse)
library(purrr)
library(fen.R.util)
library(megadaph.mtdna)

# Split the table by species, annotate each seperately, then recombine.
main <- function(var_table_file, output_file, snpeff_config) {
  var_table <- read_csv(var_table_file)
  species_tables <- split_table(var_table, by = "species")
  annotated_table <- species_tables %>%
    map(annotate_variant_table, snpeff_config = snpeff_config) %>%
    bind_rows
  write_csv(annotated_table, output_file)
}

if (!interactive()) {
  opts <- docopt::docopt(doc)
  main(
    unlist(opts["VAR_TABLE"]),
    unlist(opts["output"]),
    unlist(opts["config"]))
}
