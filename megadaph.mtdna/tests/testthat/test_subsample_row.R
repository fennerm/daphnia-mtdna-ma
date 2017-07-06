library(megadaph.mtdna)
context("Subsample row")

s <- readRDS("tests/dat/stranded_pile_row.Rds")
us <- readRDS("tests/dat/unstranded_pile_row.Rds")

s_sub <- subsample_row(s, 100)
us_sub <- subsample_row(us, 100)
