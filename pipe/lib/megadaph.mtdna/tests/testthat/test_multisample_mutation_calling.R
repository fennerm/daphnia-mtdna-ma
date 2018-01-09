library(megadaph.mtdna)
context("Multisample mutation calling")

rotation_script <- "../../exec/rotate_ref.py"
og_fa <- unlist(lapply(c(1, 2, 3, 4), function(x) {
                           file.path("..", "dat",
                                     paste("rotation_seq", x, ".fa", sep = ""))
}))

og_seqs <- lapply(og_fa, function(x) {
                    chars <- as.character(scan(x, what = character())[2])
                    as.numeric(unlist(strsplit(chars, "")))


})

rotate_ref <- function(x) {
  system(paste(rotation_script, x), intern = TRUE)
}

rot_seqs <- lapply(og_fa, function(x) {
                     as.numeric(unlist(strsplit(rotate_ref(x)[2], "")))
})

check_og_to_rot_sequence <- function(og_seq, rot_seq) {
  n <- length(og_seq)
  computed_og_idx <- unlist(lapply(1:n, rot_to_og, seq_length = n))
  expect_equal(og_seq[computed_og_idx], rot_seq)
}

check_rot_to_og_sequence <- function(og_seq, rot_seq) {
  n <- length(og_seq)
  computed_rot_idx <- unlist(lapply(1:n, og_to_rot, seq_length = n))
  expect_equal(rot_seq[computed_rot_idx], og_seq)
}

test_that("og_to_rot correct output", {
            mapply(check_og_to_rot_sequence, og_seqs, rot_seqs)
            mapply(check_rot_to_og_sequence, og_seqs, rot_seqs)
})

check_compute_split_indices <- function(og_seq, rot_seq) {
  n <- length(og_seq)
  splits <- compute_split_indices(n)
  spliced <- c(splits$og, rot_seq[splits$rot])
  expect_false(any(duplicated(spliced)))
}

test_that("compute_reference_split_indices output is correct", {
            mapply(check_compute_split_indices, og_seqs, rot_seqs)
}
)

test.tab <- readRDS("../dat/test_table.Rds")

find.test.tab <- data.frame(
  species = rep("pulex", 3),
  genotype = rep("L", 3),
  isolate = rep("L", 3),
  sample = rep(1, 3),
  pos = c(1:3),
  coverage = rep(2000, 3),
  ref = rep("A", 3),
  alt = c("A", "-", "-"),
  class = c("snv", "deletion", "deletion"),
  af = rep(1, 3),
  af.diff = c(0, 1, 1),
  strand.bias = rep(0, 3),
  unique = rep(TRUE, 3),
  coverage.prop = rep(1, 3),
  p = c(1, 0, 0)
)
