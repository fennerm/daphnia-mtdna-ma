#!/usr/bin/env Rscript
"Calculate mean sequencing depth from the spliced pileups

Usage: calc_mean_seq_depth.R --output CSV PILEUP ...

Options:
  -o --output=CSV   Output file" -> doc

library(docopt)
library(purrr)
library(tidyverse)
library(tools)

#' Calculate mean sequencing depth and standard deviation for a single pileup
calc_mean_seq_depth <- function(pileup) {
  depth <- pileup %>% pmap_int(sum)
  mean_depth <- mean(depth)
  depth_sd <- sd(depth)
  list(mean = mean_depth, sd = depth_sd)
}

main <- function(output_file, pileup_files) {
  sample <- file_path_sans_ext(basename(pileup_files))
  pileups <- pileup_files %>% map(read_csv)
  output_table <- pileups %>%
    map(calc_mean_seq_depth) %>%
    bind_rows %>%
    cbind(sample, .)
  write_csv(output_table, output_file)
}

if (!interactive()) {
  opts <- docopt(doc)
  pileup_files <- unlist(opts["PILEUP"])
  output_file <- unlist(opts["output"])
  main(output_file, pileup_files)
}
