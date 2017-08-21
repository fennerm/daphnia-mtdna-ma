library(megadaph.mtdna)
context("Highest frequency")

s <- readRDS("../dat/stranded_pile_row.Rds")
us <- readRDS("../dat/unstranded_pile_row.Rds")

s_maj <- highest_freq(s)
us_maj <- highest_freq(us)

## UNFINISHED
# test_that("Output correct", {
#     expect_equal(12, length(s_sub))
#     expect_equal(6, length(us_sub))
# })
