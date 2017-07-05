# library(Rsamtools)
# library(data.table)

## Set of functions for creating and manipulating allele count pileups.

#' @export
subsample_pileup <- function(pile, coverage_cap) {
    subsampled_rows <- apply(pile, 1, function(row) {
        row <- unlist(row)
        row_cov <- sum(row)
        if (row_cov < coverage_cap) {
            row
        } else {
            sample_pop <- unlist(sapply(1:length(row), function(i) {
                rep(i, row[i])
            }))
            samp <- sample(sample_pop, size=coverage_cap)
            sapply(1:6, function(i) length(which(samp==i)))
        }
    })
    subsampled_pile <- data.table::setDT(data.frame(t(subsampled_rows)))
    data.table::setnames(subsampled_pile, colnames(pile))
    subsampled_pile
}

## Given:
##  pile - a named numeric vector or a data.table; allele counts
## Return:
##  The name of the allele with the highest abundance at that position, or a
##  vector of major alleles for each position in the table.

#' @export
major_allele <- function(pile) {
    apply(pile, 1, function(x) {
        x <- unlist(x)
        # Get column index of major allele
        idx <- which.max(x)
        # Get allele from column index
        major_allele <- names(x)[idx]
        major_allele
    })
}

## Given:
##  x - A vector
## Return:
##  The mode
##  WARNING: Returns the first value in case of ties

#' @export
mode <- function(x) {
    ux <- unique(x)
    mode <- ux[which.max(tabulate(match(x, ux)))]
    mode
}

## Given:
##  pile - List of data.table; Each list item contains allele counts at each
##         genome position for a single sample.
##  consensus - Character vector; A multisample consensus genome sequence
## Return:
##  The consensus sequence of the samples

#' @export
create_consensus <- function(pile) {
    nsamples <- length(pile)
    ## Since all the same bases are expected to be represented in each read count
    ## file, we set the number of base pairs using the first file.
    major_alleles <- lapply(pile, major_allele)
    major_alleles <- matrix(unlist(major_alleles), ncol=nsamples)
    consensus <- apply(major_alleles, 1, mode)
    consensus
}

#' @export
create_consensus_fasta <- function(pile, filename) {
    nsamples <- length(pile)
    major_alleles <- lapply(pile, major_allele)
    major_alleles <- matrix(unlist(major_alleles), ncol = nsamples)
    consensus <- apply(major_alleles, 1, mode)
    consensus <- consensus[which(consensus != "-")]
    consensus_string <- paste0(consensus, collapse = "")
    fileconn <- file(filename)
    writeLines(consensus_string, fileconn)
    close(fileconn)
}

#' @export
piles_to_mut_wt <- function(piles, mut_consensus) {
    bp <- length(mut_consensus)
    mut_wt_counts <- lapply(piles, convert_to_mut_wt_counts, mut_consensus)
    mut_wt_matrices <- lapply(1:bp, function(i) {
        mat <- sapply(mut_wt_counts, function(x) {
            x[i,]
        })
        t(mat)
    })
    mut_wt_matrices
}

#' @export
destrand <- function(pile) {
    idx1 <- seq(1, ncol(pile)-1, by=2)
    idx2 <- idx1+1
    destranded <- mapply(function(x1, x2) {
        pile[, x1, with=FALSE] + pile[, x2, with=FALSE]
    }, idx1, idx2)
    colnames(destranded) <- c("A", "C", "G", "T", "-", "+")
    destranded
}

## Given:
##  pile - A data.table; allele counts at each genome position
##  consensus - Character vector; A multisample consensus genome sequence
## Return:
##  A two column data.frame with minor allele and its frequency.
##  Not vectorized.

#' @export
create_maf_table <- function(pile, consensus) {
    bp <- nrow(pile)

    minor_allele <- lapply(c(1:bp), function(i) {
        # Minor allele is computed by finding the major allele, excluding the
        # consensus allele.
        allele <- major_allele(pile[i, -c(consensus[i]), with=FALSE])

        # Compute allele frequency
        coverage <- sum(pile[i,])
        count <- pile[i, allele, with=FALSE]
        freq <- count/coverage
        c(allele, freq)
    })

    # Convert to data.frame and format
    minor_allele <- matrix(unlist(minor_allele), byrow=TRUE, ncol=2)
    minor_allele <- data.frame(minor_allele, stringsAsFactors=FALSE)
    colnames(minor_allele) <- c("allele", "frequency")
    minor_allele$frequency <- as.numeric(minor_allele$frequency)
    minor_allele
}

## Given:
##  pile - List of data.tables; Each list item contains allele counts at each
##         genome position for a single sample.
##  consensus - Character vector; A multisample consensus genome sequence
## Return:
##  A character vector of mutant alleles at every genome position.

#' @export
create_mut_consensus <- function(pile, consensus) {
    bp <- nrow(pile[[1]])

    # Create table of minor allele frequencies for each sample.
    maf_tables <- lapply(pile, create_maf_table, consensus)

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

## Given:
##  pile - A data.table; allele counts at each genome position
## Return:
##  A numeric vector with coverage at every position

#' @export
get_coverage <- function(pile) {

    # Function for getting coverage of a single sample.
    get_sample_coverage <- function(p) {
        unlist(apply(p, 1, sum))
    }

    # If list of data.tables, apply to each.
    if (is.null(ncol(pile))) {
        lapply(pile, get_sample_coverage)
    } else {
        get_sample_coverage(pile)
    }
}

## Given:
##  pile - A data.table; allele counts at each genome position
##  mut_consensus - Character vector; The mutant allele at each genome position
## Return:
##  A numeric vector with mutant allele counts at every position.

#' @export
get_mut_counts <- function(pile, mut_consensus) {
    bp <- nrow(pile)
    mut_count <- sapply(c(1:bp), function(i) {
        pile[i, mut_consensus[i], with=FALSE]
    })
    unlist(mut_count)
}
## Given:
##  pile - A data.table; allele counts at each genome position
##  mut_consensus - Character vector; The mutant allele at each genome position
## Return:
##  A numeric vector with mutant allele frequencies at every position.

#' @export
get_mut_af <- function(pile, mut_consensus) {
    bp <- nrow(pile)

    # Compute coverage
    coverage <- get_coverage(pile)
    # Compute frequency from counts
    mut_counts <- get_mut_counts(pile, mut_consensus)
    mut_af <- mut_counts/coverage
    mut_af
}

#' @export
stranded_indices <- function(idx) {
    str_idx <- list(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 10), c(11, 12))
    matrix(unlist(str_idx[idx]), ncol=2, byrow=TRUE)
}

#' @export
get_minor_index <- function(pile) {
    idx <- apply(pile, 1, function(row) {
        sort.int(row, index.return = TRUE, decreasing=TRUE)$ix[2]
    })
    idx
}
## Given:
##  pile - A data.table or numeric vector; allele counts
## Return:
##  Minor allele counts at each position in pile.

#' @export
get_minor_counts <- function(pile, stranded=FALSE) {
    if (stranded) {
        bp <- nrow(pile)
        ds <- destrand(pile)
        minor_idx <- get_minor_index(ds)
        stranded_idx <- stranded_indices(minor_idx)
        counts <- lapply(1:bp, function(i) {
            pile[i, stranded_idx[[i]]]
        })
    } else {
        counts <- apply(pile, 1, function(row) {
            sort.int(row)[length(row) - 1]
        })
    }
    counts
}

## Given:
##  pile - A data.table; allele counts at each genome position
##  mut_consensus - Character vector; The mutant allele at each genome position
## Return:
##  An nx2 matrix, where cols are mutant, and wild-type allele counts, and rows
##  are genome positions.

#' @export
convert_to_mut_wt_counts <- function(pile, mut_consensus) {
    mut_counts <- get_mut_counts(pile, mut_consensus)
    covs <- get_coverage(pile)
    non_mut_counts <- covs - mut_counts
    cbind(mut_counts, non_mut_counts)
}

## Given:
##  pile - A data.table; allele counts at each genome position
##  mut_consensus - Character vector; The mutant allele at each genome position
## Return:
##  An nx2 matrix, where cols are mutant allele counts, and coverage, and rows
##  are genome positions.

#' @export
convert_to_mut_cov_counts <- function(pile, mut_consensus) {
    mut_counts <- get_mut_counts(pile, mut_consensus)
    depths <- get_coverage(pile)
    cbind(mut_counts, depths)
}

#' @export
create_pileup <- function(bam, min_base_quality = 30,
                          distinguish_strands = FALSE) {
    pileup_param <- Rsamtools::PileupParam(
        max_depth = 1000000,
        distinguish_strands = distinguish_strands,
        distinguish_nucleotides = TRUE,
        ignore_query_Ns = TRUE,
        min_nucleotide_depth = 0,
        include_deletions = TRUE,
        include_insertions = TRUE,
        min_base_quality = min_base_quality)

    pile <- data.table::setDT(Rsamtools::pileup(bam,
                                                pileupParam = pileup_param))

    if (distinguish_strands) {
        pile_wide <- data.table::dcast(pile,
                                       seqnames+pos~nucleotide+strand,
                                       value.var = "count")
    } else {
        pile_wide <- data.table::dcast(pile,
                                          seqnames+pos ~ nucleotide,
                                          value.var = "count")
    }
    pile_wide <- pile_wide[, 3:length(pile_wide)]

    for (j in seq_len(ncol(pile_wide)))
        set(pile_wide,which(is.na(pile_wide[[j]])),j,0)
    pile_wide
}


## Given:
## og, rot - character vector; filepaths to OG and ROT BAM files
## min_base_quality - numeric; minimum Phred base quality for base to be
##                    included in table
#' @export
create_spliced_count_table <- function(og, rot, min_base_quality=30,
                                       distinguish_strands=FALSE) {
    species <- get_species(og)
    # Parameters for creating pileup from bam files.
    og_pile <- create_pileup(og, min_base_quality, distinguish_strands)
    rot_pile <- create_pileup(rot, min_base_quality, distinguish_strands)

    # Splice OG and ROT tables
    pile_wide <- splice_mtdna_data(og_pile, rot_pile, species)

    pile_wide
}

