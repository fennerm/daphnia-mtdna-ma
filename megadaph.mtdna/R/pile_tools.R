#' Set of functions for creating and manipulating allele count pileups.


#' Subsample a count vector
#' @param r An integer vector
#' @param cap Integer; Vector will be subsampled such that its sum = cap
#' @return An integer vector
#' @export
subsample <- function(r, cap) {
    r <- unlist(r)
    cov <- sum(r)
    if (cov > cap) {
        # Calculate the surplus and randomly subtract from the row
        surplus <- cov - cap
        prop <- r / sum(r)
        smp <- rmultinom(1, surplus, prop)
        new_row <- as.integer(r - smp)
    } else {
        new_row <- r
    }
    new_row
}

#' Subsample an allele count pileup
#' @param pile A pileup
#' @param cap Integer; Pileup will be subsampled such that its sum = cap
#' @return A pileup
#' @export
#' @importFrom data.table setDT setnames
subsample_pileup <- function(pile, cap) {
    subsampled_rows <- apply(pile, 1, subsample)
    subsampled_pile <- setDT(data.frame(t(subsampled_rows)))
    setnames(subsampled_pile, colnames(pile))
    subsampled_pile
}

#' Calc. the mode of a vector.
#'
#' Returns the first value in case of ties
#' @export
mode <- function(x) {
    ux <- unique(x)
    mode <- ux[which.max(tabulate(match(x, ux)))]
    mode
}

#' Create a consensus sequence from a set of nucleotide count pileups
#' @param piles List of data.table; Each list item contains allele counts at
#'        each genome position for a single sample.
#' @param fasta File path as character; If given, consensus sequence will be
#'        written to the path
#' @return Character vector; The multisample consensus sequence as a string
#' @export
create_consensus <- function(piles, fasta=NULL) {
    nsamples <- length(piles)
    #' Determine the highest frequency allele at each site for each sample
    major_alleles <- lapply(piles, major_allele)
    #' Since all the same bases are expected to be represented in each read
    #' count file, we set the number of base pairs using the first file.
    major_alleles <- matrix(unlist(major_alleles), ncol = nsamples)
    consensus <- apply(major_alleles, 1, mode)

    if (!is.null(fasta)) {
        consensus_string <- paste0(consensus, collapse = "")
        fileconn <- file(fasta)
        writeLines(consensus_string, fileconn)
        close(fileconn)
    }

    consensus
}


#' Splice the original and rotated BAM files into a single allele count table
#' @param og Character vector; filepath to the original (og) .bam file
#' @param rot Character vector; filepath to the rotated (rot) .bam file
#' @param min_base_quality Numeric; minimum Phred base quality for base to be
#'        included in table
#' @param distinguish_strands; If TRUE base counts on the +/- strands are
#'        counted separately
#' @return A data.frame
#' @export
create_spliced_count_table <- function(og, rot, min_base_quality=30,
                                       distinguish_strands=FALSE) {
    species <- get_species(og)
    # Parameters for creating pileup from bam files.
    og_pile <- create_pileup(og, min_base_quality, distinguish_strands)
    rot_pile <- create_pileup(rot, min_base_quality, distinguish_strands)

    # Splice OG and ROT tables
    spliced_pile <- splice_mtdna_data(og_pile, rot_pile, species)
    spliced_pile
}

#' Create a sequence of the consensus alternate (mutant) alleles at each genome
#' position among multiple samples.
#'
#' The alternate allele is the non-consensus base with the highest frequency at
#' a locus within a single sample. These alternate alleles may be true mutants
#' or might just be sequencing error.
#' @param piles List of data.table; Each list item contains allele counts at
#'        each genome position for a single sample.
#' @param consensus character; The consensus sequence for the samples
#' @return character; The alternate consensus sequence
#' @export
create_mut_consensus <- function(piles, consensus) {
    bp <- nrow(piles[[1]])

    # Create table of minor allele frequencies for each sample.
    maf_tables <- lapply(piles, create_maf_table, consensus)

    # Create vector with index of the sample with the largest minor allele
    # frequency at each position.
    max_idx <- sapply(c(1:bp), function(i) {
        which.max(sapply(maf_tables, function(x) {
            x[i, 2]
        }))
    })

    # Create a vector with minor alleles of samples in max_idx
    mut_consensus <- sapply(c(1:bp), function(i) {
        maf_tables[[max_idx[i]]][i, 1]
    })
    mut_consensus
}

#' Convert allele count pileups to mutant vs wildtype pileups
#'
#' @param piles List of data.table; Each list item contains allele counts at
#'        each genome position for a single sample.
#' @param mut_consensus; An alternate allele consensus sequence produced by
#'        `create_mut_consensus`
#' @export
piles_to_mut_wt <- function(piles, mut_consensus) {
    bp <- length(mut_consensus)
    mut_wt_counts <- lapply(piles, convert_to_mut_wt_counts, mut_consensus)
    mut_wt_matrices <- lapply(1:bp, function(i) {
        mat <- sapply(mut_wt_counts, function(x) {
            x[i, ]
        })
        t(mat)
    })
    mut_wt_matrices
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

#' Get the highest frequency allele in a vector of allele counts
#' @export
highest_freq <- function(v) {
    n <- length(v)
    v <- unlist(v)
    if (n == 12) {
        v <- destrand(v)
    }

    # Get column index of major allele
    idx <- which.max(v)
    # Get allele from column index
    major_allele <- names(v)[idx]
    major_allele
}

#' Get the major allele at each genome position for a single sample
#'
#' WARNING: Returns the first value if tied
#' @param pile - A data.table; allele counts at each genome position
#' @return Character; Vector of major alleles
#' @export
major_allele <- function(pile) {
    apply(pile, 1, highest_freq)
}

#' Generate a table of minor allele frequencies
#' @param pile - A data.table; allele counts at each genome position
#' @param consensus character; The consensus sequence for the samples
#' @return data.table; Two columns - minor allele and minor allele frequency
#' @export
create_maf_table <- function(pile, consensus) {
    bp <- nrow(pile)

    minor_alleles <- lapply(c(1:bp), function(i) {
        # Minor allele is computed by finding the major allele, excluding the
        # consensus allele.
        allele <- major_allele(pile[i, -c(consensus[i]), with = FALSE])

        # Compute allele frequency
        coverage <- sum(pile[i, ])
        count <- pile[i, allele, with = FALSE]
        freq <- count / coverage
        c(allele, freq)
    })

    # Convert to data.frame and format
    minor_alleles <- matrix(unlist(minor_alleles), byrow = TRUE, ncol = 2)
    minor_alleles <- data.frame(minor_alleles, stringsAsFactors = FALSE)
    colnames(minor_alleles) <- c("allele", "frequency")
    minor_alleles$frequency <- as.numeric(minor_alleles$frequency)
    minor_alleles
}


#' Get sequencing coverage for each genome position
#' @param pile - A data.table; allele counts at each genome position
#' @return A numeric vector
#' @export
get_coverage <- function(pile) {
    # Function for getting coverage of a single sample.
    get_sample_coverage <- function(p) {
        unlist(apply(pile, 1, sum))
    }

    # If list of data.tables, apply to each.
    if (is.null(ncol(pile))) {
        lapply(pile, get_sample_coverage)
    } else {
        get_sample_coverage(pile)
    }
}

#' Get mutant allele counts at each genome position
#' @param pile - A data.table; allele counts at each genome position
#' @param mut_consensus; An alternate allele consensus sequence produced by
#'        `create_mut_consensus`
#' @return Numeric vector
#' @export
get_mut_counts <- function(pile, mut_consensus) {
    mut_counts <- sapply(c(1:nrow(pile)), function(i) {
        pile[i, mut_consensus[i], with = FALSE]
    })
    unlist(mut_counts)
}

#' Get mutant allele frequency at each genome position
#' @param pile - A data.table; allele counts at each genome position
#' @param mut_consensus; An alternate allele consensus sequence produced by
#'        `create_mut_consensus`
#' @param Numeric vector
#' @export
get_mut_af <- function(pile, mut_consensus) {
    # Compute coverage
    coverage <- get_coverage(pile)
    # Compute frequency from counts
    mut_counts <- get_mut_counts(pile, mut_consensus)
    mut_af <- mut_counts / coverage
    mut_af
}

#' Given a column index into an unstranded pileup, return the equiv. stranded
#' indices
.stranded_indices <- function(idx) {
    str_idx <- list(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 10), c(11, 12))
    str_idx <- matrix(unlist(str_idx[idx]), ncol = 2, byrow = TRUE)
    str_idx
}

#' Convert an allele count pileup to a table of mutant vs. wt counts
#' @param pile - A data.table; allele counts at each genome position
#' @param mut_consensus; An alternate allele consensus sequence produced by
#'        `create_mut_consensus`
#' @return An nx2 matrix, where cols are mutant, and wild-type allele counts,
#' and rows are genome positions.
#' @export
convert_to_mut_wt_counts <- function(pile, mut_consensus) {
    mut_counts <- get_mut_counts(pile, mut_consensus)
    covs <- get_coverage(pile)
    non_mut_counts <- covs - mut_counts
    cbind(mut_counts, non_mut_counts)
}

#' Convert an allele count pileup to a table of mutant vs. coverage
#' @param pile - A data.table; allele counts at each genome position
#' @param mut_consensus; An alternate allele consensus sequence produced by
#'        `create_mut_consensus`
#' @return An nx2 matrix, where cols are mutant allele counts, and coverage, and
#' rows are genome positions.
#' @export
convert_to_mut_cov_counts <- function(pile, mut_consensus) {
    mut_counts <- get_mut_counts(pile, mut_consensus)
    depths <- get_coverage(pile)
    cbind(mut_counts, depths)
}

#' Create an allele count pileup with RSamtools
#' @param bam Path to a .bam file
#' @param min_base_quality All read positions with phred score <
#'        min_base_quality will be excluded
#' @param distinguish_strands; If TRUE, base counts on the +/- strands are
#'        counted separately
#' @param min_nucleotide_depth; Minimum allele count to be included in pileup
#' @importFrom Rsamtools PileupParam pileup
#' @importFrom data.table setDT
#' @importFrom reshape2 dcast
#' @export
create_pileup <- function(bam, min_base_quality = 30,
                          distinguish_strands = FALSE,
                          min_nucleotide_depth = 1) {
    pileup_param <- PileupParam(
        max_depth = 1000000,
        distinguish_strands = distinguish_strands,
        distinguish_nucleotides = TRUE,
        ignore_query_Ns = TRUE,
        min_nucleotide_depth = min_nucleotide_depth,
        include_deletions = TRUE,
        include_insertions = TRUE,
        min_base_quality = min_base_quality)

    cat("Rsamtools pileup \n")
    pile <- pileup(bam, pileupParam = pileup_param)

    #' Data.table dcast has a known bug which causes errors if pile is large.
    cat("Casting to wide format \n")
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
    cat("Removing extra cols \n")
    pile_wide[, c("seqnames", "pos") := NULL]

    pile_wide
}
