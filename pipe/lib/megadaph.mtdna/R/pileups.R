#' Splice the original and rotated BAM files into a single allele count table
#' @param og_bam .bam files from alignment to the unrotated reference
#'                sequences
#' @param rot_bam .bam files from alignment to the rotated reference
#'                sequences
#' @param min_base_quality Numeric; minimum Phred base quality for base to be
#'        included in table
#' @param distinguish_strands; If TRUE base counts on the +/- strands are
#'        counted separately
#' @return A data.frame
#' @export
construct_spliced_pileup <- function(og_bam, rot_bam, min_base_quality=30,
                                     distinguish_strands=FALSE) {
  # Create individual pileups
  og_pile <- construct_pileup(og_bam, min_base_quality, distinguish_strands)
  rot_pile <- construct_pileup(rot_bam, min_base_quality, distinguish_strands)

  # Splice the pileups
  spliced_pile <- splice_pileups(og_pile, rot_pile, species)
  spliced_pile
}

#' Create an allele count pileup with RSamtools
#' @param bam Path to a .bam file
#' @param min_base_quality All read positions with phred score <
#'        min_base_quality will be excluded
#' @param distinguish_strands; If TRUE, base counts on the +/- strands are
#'        counted separately
#' @param min_nucleotide_depth; Minimum allele count to be included in pileup
#' @return A data.table
#' @importFrom Rsamtools PileupParam pileup
#' @importFrom data.table setDT
#' @importFrom reshape2 dcast
#' @export
construct_pileup <- function(bam, min_base_quality = 30,
                             distinguish_strands = FALSE,
                             min_nucleotide_depth = 1) {
  pileup_param <- PileupParam(max_depth = 1000000,
                              distinguish_strands = distinguish_strands,
                              distinguish_nucleotides = TRUE,
                              ignore_query_Ns = TRUE,
                              min_nucleotide_depth = min_nucleotide_depth,
                              include_deletions = TRUE,
                              include_insertions = TRUE,
                              min_base_quality = min_base_quality)

  pile <- pileup(bam, pileupParam = pileup_param)

  # Cast the table to wide format
  if (distinguish_strands) {
    pile_wide <- dcast(pile,
                       seqnames + pos ~ nucleotide + strand,
                       value.var = "count",
                       fill = 0L)
  } else {
    pile_wide <- dcast(pile,
                       seqnames + pos ~ nucleotide,
                       value.var = "count",
                       fill = 0L)
  }
  pile_wide <- setDT(pile_wide)
  # Remove unnecessary columns
  pile_wide[, c("seqnames", "pos") := NULL]
  pile_wide
}

#' Consolidate strand counts in a strand-split allele count table
#' @param counts data.table; Table with nucleotide and indel counts for each
#'        genome position split by strand
#' @return data.table; A allele count table with the stranded counts summed
#' @export
#' @importFrom data.table setDT
destrand <- function(counts) {
  idx1 <- seq(1, 11, by = 2)
  idx2 <- idx1 + 1

  if (is.numeric(counts)) {
    destranded <- counts[idx1] + counts[idx2]
  } else {
    sum_counts <- function(x1, x2) {
      counts[, x1, with = FALSE] + counts[, x2, with = FALSE]
    }
    destranded <- mapply(sum_counts, idx1, idx2)
    destranded <- setDT(destranded)
  }
  names(destranded) <- c("A", "C", "G", "T", "-", "+")
  destranded
}

#' Splice the 'original' and 'rotated' pileups into a single pileup
#' @param og_pile Pileup from alignment to the 'original' reference sequence
#' @param rot_pile Pileup from alignment to the 'rotated' reference sequence
#' @return A spliced pileup containing the middle halves of both the 'rotated'
#' and 'original' pileups
#' @importFrom data.table rbindlist
splice_pileups <- function(og_pile, rot_pile, rot_ref) {
  bp <- nrow(og_pile)
  splx <- compute_split_indices(bp)

  data_type <- class(og_pile)[1]

  if (data_type == "matrix") {
    spliced <- rbind(rot_pile[splx[[1]], ], og_pile[splx[[2]], ],
                     rot_pile[splx[[3]], ])
  } else if (data_type == "data.table") {
    spliced <- rbindlist(list(rot_pile[splx[[1]], ], og_pile[splx[[2]], ],
                              rot_pile[splx[[3]], ]))
  } else if (data_type %in% c("numeric", "integer")) {
    spliced <- c(rot_pile[splx[[1]]], og_pile[splx[[2]]], rot_pile[splx[[3]]])
  }
  spliced
}

#' Compute the indices for splicing the 'original' and 'rotated' pileups
#' @param sequence_length Integer; The length of the reference sequence
#' @return List, length=3; First and third list items are indices from 'rotated'
#' pileup. Second item is indices into 'original' pileup.
compute_split_indices <- function(sequence_length) {
  spl <- floor(seq(1, sequence_length, len = 5))
  rot1 <- og_to_rot(spl[1], sequence_length):og_to_rot(spl[2] - 1,
                                                       sequence_length)
  og <- spl[2]:(spl[4] - 1)
  rot2 <- og_to_rot(spl[4], sequence_length):og_to_rot(spl[5], sequence_length)
  list(rot1, og, rot2)
}

