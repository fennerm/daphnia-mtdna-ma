#!/usr/bin/env Rscript
## Create consensus sequences for Sanger sequencing of mtDNA control region

INPUT_DIR <- file.path('..', '..', 'data', 'sanger_dloop_nov_2017')
AB1_DIR <- file.path(INPUT_DIR, 'ab1')
CHROMATOGRAMS <- list.files(AB1_DIR, full.names = TRUE)
NAMES <- c('FBSC-Control Region-Large Sanger Product',
           'GASC-Control Region-Small Sanger Product',
           'GASC-Control Region- Large Sanger Product')
OUTPUT_DIR <- file.path(INPUT_DIR, 'merged_seqs')

#' Pair the chromatogram by sample
#' @param chromatograms A list of chromatograms
#' @return A named list of lists. First list is forward sequences,
#'         second is reverse.
.pair_chromatograms <- function(chromatograms) {
    split_list <- split(chromatograms, 1:2)
    names(split_list) <- c('fwd', 'rev')
    split_list
}


#' Create a single consensus sequence from forward and reverse sequences
#' @param fwd,rev Forward and reverse chromatograms
#' @importFrom Biostrings DNAStringSet reverseComplement
#' @importFrom sangeranalyseR merge.reads
#' @importFrom sangerseqR primarySeq readsangerseq
#' @export
create_consensus <- function(fwd, rev) {
    # Get sequence
    fwd <- primarySeq(readsangerseq(fwd))
    rev <- primarySeq(readsangerseq(rev))
    rev <- reverseComplement(rev)

    # this gives us an unaligned set of the reads we wish to merge
    reads <- DNAStringSet(c(as.character(fwd), as.character(rev)))
    names(reads) <- c('fwd', 'rev')
    merged_reads <- merge.reads(reads)
    merged_reads
}


#' Write the consensus to fasta file
#' @param sequence The fasta sequence
#' @param name The name to write to the ID field
#' @param filename The path to the output fasta
#' @export
write_fasta_sequence <- function(sequence, name, filename) {
    id <- paste('>', name)
    cat(c(id, sequence), file = filename, append = TRUE, sep = '\n')
}

#' Main function
#' @param fwd,rev Forward and reverse chromatograms
#' @param name Name of the sample (for writing fasta ID)
main <- function(fwd, rev, name) {
    merged <- create_consensus(fwd, rev)
    consensus <- as.character(merged$consensus)
    dir.create(OUTPUT_DIR, showWarnings = FALSE)
    filename <- file.path(OUTPUT_DIR, paste(name, '.fasta', sep = ''))
    if (file.exists(filename)) file.remove(filename)
    write_fasta_sequence(consensus, name, filename)
}

if (!interactive()) {
    paired <- .pair_chromatograms(CHROMATOGRAMS)
    merged <- mapply(main, fwd = paired$fwd, rev = paired$rev, name = NAMES, SIMPLIFY = FALSE)
}
