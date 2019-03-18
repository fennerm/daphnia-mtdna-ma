BASES <- c("A", "C", "G", "T")

#' Generate an empty square matrix
#'
#' @param dimnames character; The column/row dimensions.
#' @return A matrix with dimnames as the column and row names.
#'
#' @export
gen_square_matrix <- function(dimnames, fill = NA) {
  dim = length(dimnames)
  matrix(data = fill, nr = dim, nc = dim, dimnames = rep(list(dimnames), 2))
}


#' Count the number of transitions in a flat list of snv counts
#'
#' @export
count_transitions <- function(counts) {
  counts["A→G"] + counts["G→A"] + counts["C→T"] + counts["T→C"]
}


#' Count the number of transversions in a flat list of snv counts
#'
#' @export
count_transversions <- function(counts) {
  counts["A→C"] + counts["C→A"] + counts["A→T"] + counts["T→A"] +
    counts["G→C"] + counts["C→G"] + counts["G→T"] + counts["T→G"]
}

#' Generate a 4x4 matrix containing counts of base substitution types
#'
#' @param flat logical; If TRUE, return a named numeric vector instead
#' @return Matrix or list
#' @importFrom gtools permutations
#' @import tidyverse
#' @export
extract_base_sub_spectrum <- function(variant_table, flat = FALSE) {
  spectrum <- gen_square_matrix(BASES, fill = 0)
  ref_bases <- unlist(variant_table[, "ref"])
  alt_bases <- unlist(variant_table[, "alt"])
  for (i in 1:length(ref_bases)) {
    if (!is.null(ref_bases[i])) {
      if (ref_bases[i] %in% BASES && alt_bases[i] %in% BASES) {
        spectrum[ref_bases[i], alt_bases[i]] <- (
          spectrum[ref_bases[i], alt_bases[i]] + 1
        )
      }
    }
  }
  if (flat == TRUE) {
    base_perm <- permutations(4, 2, BASES)
    sub_names <- paste0(base_perm[, 1], "→", base_perm[, 2])
    spectrum <- map2_dbl(base_perm[, 1], base_perm[, 2], function(i, j) {
      spectrum[i, j]
      }
    )
    names(spectrum) <- sub_names
  }
  spectrum
}


#' Count number of items in `positions` within `interval`
num_within_interval <- function(interval, positions) {
  length(which(positions %in% interval))
}

# #'
# #' @param variant_table
# #' @importFrom fen.R.util split_into_chunks
# #' @return
# uniform_distribution_test <- function(
#   variant_table,
#   reference_length,
#   num_intervals
# ) {
#   intervals <- split_into_chunks(1:reference_length, num_intervals)
#   interval_counts <- intervals %>%
#     map(~num_within_interval(., unlist(variant_table$pos))) %>%
#     unlist
#   chisq.test(interval_counts, names(interval_counts))
# }
