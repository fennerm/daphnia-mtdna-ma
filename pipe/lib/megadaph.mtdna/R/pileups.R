#' Splice the original and rotated BAM files into a single allele count table
#' @param og_bam .bam files from alignment to the unrotated reference
#'                sequences
#' @param rot_bam .bam files from alignment to the rotated reference
#'                sequences
#' @param bp Number of base pairs in the mitochondrial reference sequence.
#' @param min_base_quality Numeric; minimum Phred base quality for base to be
#'        included in table
#' @param distinguish_strands; If TRUE base counts on the +/- strands are
#'        counted separately
#' @return A matrix
#' @export
construct_spliced_pileup <- function(og_bam, rot_bam, bp, min_base_quality = 30,
                                     distinguish_strands = FALSE) {
  # Create individual pileups
  og_pile <- construct_pileup(og_bam, bp, min_base_quality, distinguish_strands)
  rot_pile <- construct_pileup(rot_bam, bp, min_base_quality,
                               distinguish_strands)

  # Splice the pileups
  spliced_pile <- splice_pileups(og_pile, rot_pile, species)
  spliced_pile
}

#' Create an allele count pileup with RSamtools
#' @param bam Path to a .bam file
#' @param bp Number of base pairs in the mitochondrial reference sequence.
#' @param min_base_quality All read positions with phred score <
#'        min_base_quality will be excluded
#' @param distinguish_strands; If TRUE, base counts on the +/- strands are
#'        counted separately
#' @return A matrix
#' @importFrom Rsamtools PileupParam pileup
#' @importFrom reshape2 dcast
#' @importFrom fen.R.util add_missing_columns
#' @export
construct_pileup <- function(bam, bp, min_base_quality = 30,
                             distinguish_strands = FALSE) {
  pileup_param <- PileupParam(max_depth = 1000000,
                              distinguish_strands = distinguish_strands,
                              distinguish_nucleotides = TRUE,
                              min_nucleotide_depth = 1,
                              include_deletions = TRUE,
                              include_insertions = TRUE,
                              min_base_quality = min_base_quality)

  pile <- pileup(bam, pileupParam = pileup_param)
  seqname <- unique(pile$seqnames)

  # Add missing rows to the pileup
  tmp_pile <- data.frame(pos = setdiff(1:bp, pile$pos))
  # Remove NA values
  pile <- merge(pile, tmp_pile, all = TRUE)
  pile$seqnames[is.na(pile$seqnames)] <- seqname
  pile$strand[is.na(pile$strand)] <- "+"
  pile$nucleotide[is.na(pile$nucleotide)] <- "A"


  # Cast the table to wide format
  if (distinguish_strands) {
    pile$count[is.na(pile$count)] <- 0
    pile_wide <- dcast(pile,
                       seqnames + pos ~ nucleotide + strand,
                       value.var = "count",
                       fill = 0L)
    expected_colnames <- c("A_+", "A_-", "C_+", "C_-", "G_+", "G_-", "T_+",
                           "T_-", "+_-", "+_+", "-_-", "-_+")
  } else {
    pile_wide <- dcast(pile,
                       seqnames + pos ~ nucleotide,
                       value.var = "count",
                       fill = 0L)
    expected_colnames <- c("A", "C", "G", "T", "+", "-")
  }
  # Remove unnecessary columns
  pile_wide[, c("seqnames", "pos")] <- NULL

  pile_wide <- as.matrix(pile_wide)

  # Ensure that the output pileup has all expected column names
  pile_wide <- add_missing_columns(pile_wide, expected_colnames, fill = 0)
  pile_wide
}

#' Consolidate strand counts in a strand-split allele count table
#' @param pile matrix; Table with nucleotide and indel pile for each
#'        genome position split by strand
#' @return Matrix; A allele count table with the stranded pile summed
#' @export
destrand <- function(pile) {
  idx1 <- seq(1, 11, by = 2)
  idx2 <- idx1 + 1

  destranded <- pile[, idx1] + pile[, idx2]
  colnames(destranded) <- c("A", "C", "G", "T", "-", "+")
  destranded
}

#' Splice the 'original' and 'rotated' pileups into a single pileup
#' @param og_pile Pileup from alignment to the 'original' reference sequence
#' @param rot_pile Pileup from alignment to the 'rotated' reference sequence
#' @return A spliced pileup containing the middle halves of both the 'rotated'
#' and 'original' pileups
splice_pileups <- function(og_pile, rot_pile, rot_ref) {
  bp <- nrow(og_pile)
  splx <- compute_split_indices(bp)
  spliced <- rbind(rot_pile[splx[[1]], ], og_pile[splx[[2]], ],
                   rot_pile[splx[[3]], ])
  spliced
}

#' Compute the indices for splicing the 'original' and 'rotated' pileups
#'
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


#' Downsample a pileup to an even and uniform sequencing depth.
#'
#' @export 
downsample_pileup <- function(pileup, depth) {
  downsampled <- t(apply(pileup, 1, function(x) downsample_vector(x, depth)))
  colnames(downsampled) <- colnames(pileup)
  downsampled
}


#' Downsample a vector of allele counts to a target sequencing depth.
#'
#' @importFrom purrr map
downsample_vector <- function(vec, depth) {
  if (sum(vec) > depth) {
    allele_pop <- unlist(map(list(vec[vec > 0]), ~rep(names(.), .)))
    downsampled_pop <- sample(allele_pop, depth)
    downsampled <- map_int(names(vec), ~sum(downsampled_pop == .))
  } else {
    downsampled <- vec
  }
  downsampled
}
