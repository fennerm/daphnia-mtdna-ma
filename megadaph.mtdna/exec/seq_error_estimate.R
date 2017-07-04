#!/usr/bin/env Rscript

## Script for calculating sequencing error rates for heterozygous samples.
## First argument is directory location containing bam files resultant from
## mapping to nuclear genome.
## Second argument is the output file location (output is a .csv file)

# Calculate the sequencing error rate for a single bam file.

#' @export
calc_seq_err <- function(path, pileup_param) {
    # Create base count pilep and convert to a data.table
    pile <- setDT(pileup(path, pileupParam = pileup_param))

    # Convert from long to wide format and exclude position and c'some fields
    pile_wide <- dcast(
        pile, seqnames+pos ~ nucleotide, value.var="count")
    pile_wide <- pile_wide[, 3:length(pile_wide)]

    # Remove NAs
    for (j in seq_len(ncol(pile_wide)))
        set(pile_wide,which(is.na(pile_wide[[j]])),j,0)

    # Remove pile from memory
    rm(pile)

    base_pile <- pile_wide[, 1:4]
    # Count of minor allele at each position
    minor_counts <- apply(base_pile, 1, get_minor_counts)
    # Coverage at each position
    sum_counts <- apply(pile_wide, 1, sum)

    # Homozygous positions
    #0.2 chosen as cutoff since pbinom is approx 0.001 for depth of 20.
    cutoff <- 0.2 * sum_counts
    homo <- which((minor_counts < cutoff) &
                      (sum_counts > 39) &
                  (unlist(pile_wide[, "+"]) < cutoff) &
                  (unlist(pile_wide[, "-"]) < cutoff))

    # Sequencing error estimate
    seq_err <- mean(minor_counts[homo] / sum_counts[homo])
    # Insertion error estimate
    ins_err <- mean(unlist(pile_wide[homo, "+"])/sum_counts[homo])
    # Deletion error estimate
    del_err <- mean(unlist(pile_wide[homo, "-"])/sum_counts[homo])

    err <- c(seq_err, ins_err, del_err)
    err
}

# Calculate the sequencing error rate for a list of files

#' @export
seq_err <- function(file_list, pileup_param) {

    # Only use multithreading if run as script
    if(interactive()) {
        seq_errs <- lapply(file_list, calc_seq_err, pileup_param)
    } else {
	    ncores <- length(file_list)
	    cl <- parallel::makeCluster(ncores)
        parallel::clusterEvalQ(cl, library(data.table))
        parallel::clusterEvalQ(cl, library(Rsamtools))
	    parallel::clusterExport(cl, "get_minor_counts")
        seq_errs <- parallel::parLapply(cl, file_list, calc_seq_err, pileup_param)
        parallel::stopCluster(cl)
    }

    # Convert to a matrix and name cols and rows
    seq_errs <- matrix(unlist(seq_errs), ncol=3, byrow=TRUE)
    colnames(seq_errs) <- c("sequencing_error", "insertion_error",
                            "deletion_error")
    rownames(seq_errs) <- basename(file_list)
    seq_errs
}

# Parameters for creating pileup from bam files.
pileup_param <- Rsamtools::PileupParam(max_depth=1000000,
                            distinguish_strands=FALSE,
                            distinguish_nucleotides = TRUE,
                            ignore_query_Ns = TRUE,
                            include_deletions = TRUE,
                            include_insertions = TRUE)

# If run as script, use command line arguments. The script also works fine in
# interactive mode - just need to define a character vector of filenames and
# then run seq_err function on it.
# if (!interactive()) {
#     args = commandArgs(trailingOnly=TRUE)
#     file_list <- list_bams(args[1])
#     seq_err_table <- seq_err(file_list, pileup_param)
#     write.csv(seq_err_table, args[2], quote=FALSE)
# }
