library(megadaph.mtdna)
context("Subsample row")

s <- readRDS("tests/dat/stranded_pile_row.Rds")
us <- readRDS("tests/dat/unstranded_pile_row.Rds")

s_sub <- subsample_row(s, 100)
us_sub <- subsample_row(us, 100)

test_that("Output length", {
    expect_equal(12, length(s_sub))
    expect_equal(6, length(us_sub))
})

test_that("Sums to coverage cap", {
    expect_equal(sum(s_sub), 100)
    expect_equal(sum(us_sub), 100)
})

test_that("Returns a vector", {
    expect_equal(class(s_sub), "integer")
    expect_equal(class(us_sub), "integer")

})
