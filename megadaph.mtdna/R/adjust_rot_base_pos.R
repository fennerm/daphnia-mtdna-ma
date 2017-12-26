#!/usr/bin/env Rscript

#' Convert indices in rotated fasta file into the pre-rotation indices.
#' Under rotation, the second half of fasta sequence is moved to start.
#' See rotate_ref.py
#' @param pos Index in rotated fasta file
#' @param seq_length The length of the fasta sequence
#' @return The rotated index
#' @export
rot_to_og <- function(pos, seq_length) {
    if ((pos > seq_length) || (pos < 0)) {
        stop("Invalid position")
    } else {
        midpoint <- round(seq_length / 2)
        odd_even <-  seq_length %% 2
        if (pos < (midpoint + 1 - odd_even)) {
            pos + midpoint
        } else {
            pos - midpoint + odd_even
        }
    }
}

#' Convert indices in original fasta file into rotated indices.
#' Under rotation, the second half of fasta sequence is moved to start.
#' @param pos Index in original fasta file
#' @param seq_length The length of the fasta sequence
#' @return The rotated index
#' @export
og_to_rot <- function(pos, seq_length) {
    if ((pos > seq_length) || (pos < 0)) {
        stop("Invalid position")
    } else {
        midpoint = round(seq_length / 2)
        ## We need to adjust by 1 if the sequence length is odd
        odd_even <-  seq_length %% 2
        if (pos < (midpoint + 1)) {
            pos + midpoint - odd_even
        } else {
            pos - midpoint
        }
    }
}
