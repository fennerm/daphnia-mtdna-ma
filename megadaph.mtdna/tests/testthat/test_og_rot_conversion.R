library(megadaph.mtdna)
context("OG <-> ROT conversion")

## Algorithm used for rotating the sequence - Included for convenience
rotate_alg <- function(seq) {
    bp <- length(seq)
    midpoint <- round(bp / 2)
    first_half <- seq[1:(midpoint)]
    second_half <- seq[(midpoint + 1):bp]
    c(second_half, first_half)
}

test_that("Correct output", {
    expect_equal(og_to_rot(5, 10), 10)
    expect_equal(og_to_rot(6, 10), 1)
    expect_equal(og_to_rot(6, 11), 11)
    expect_equal(og_to_rot(7, 11), 1)
    expect_equal(rot_to_og(5, 10), 10)
    expect_equal(rot_to_og(6, 10), 1)
    expect_equal(rot_to_og(5, 11), 11)
    expect_equal(rot_to_og(6, 11), 1)
})
