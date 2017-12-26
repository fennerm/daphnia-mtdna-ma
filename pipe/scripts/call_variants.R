#' Variants were called using a novel contingency matrix approach.
#' For each group of samples with a common ancestor genotype:
#' 1. Determine the major and minor allele at every locus in each sample.
#'    - The major allele is defined as the nucleotide which is present in the
#'      the most reads at a site.
#' 2. Determine the consensus major allele at each site across all samples.
#' 3. Define the minor allele as the allele which reached the highest allele
#'    frequency in a single sample, other than the major allele.
#' 4. For each site in the genome, produce a contingency table of the form
#'              | Major | Minor |
#'     Sample 1 |   -   |   -   |
#'        ...   |   -   |   -   |
#'     Sample N |   -   |   -   |
#'
#'    where 'Major' and 'Minor' are the major and minor allele frequency in
#'    each sample.
#' 5. Test the null hypothesis that the allele frequencies at a particular site
#'    are equal across samples.
#'    - For this purpose we adapted a statistical test based on the expected
#'      distribution of the maximum of a hypergeometric distribution:
#'      The exact distribution of the maximum, minimum and the range of
#'      Multinomial/Dirichlet and Multivariate Hypergeometric frequencies,
#'      Corrado, 2011

#' Process Bam files from a single isolate and construct pileups,
#' mut-wt-matrices and a 'test_table' with a lot of info relevant to variant
#' calling.
#' @param og_bams List of bam files from alignment to the unrotated reference
#'                sequences
#' @param rot_bams List of bam files from alignment to the rotated reference
#'                sequences
#' @param seq_err_rates List of sequencing error rate estimates for each sample
#' @return A three element list: [1]. A list of the pileups for each sample in
#'         isolate. [2]. A list of mut-wt contingency matrices for each sample
#'         [3]. A single data.frame 'test_table' with variant info at every
#'         genome position.
#' @export
construct_tables <- function(og_bams, rot_bams, seq_err) {
    og_bams <- unlist(og_bams)
    rot_bams <- unlist(rot_bams)
    seq_err <- unlist(seq_err)

    # Number of samples for isolate
    nsamples <- length(og_bams)

    # Create allele count pileups spliced across the OG and ROT alignments for
    # each sample.
    piles <- lapply(1:nsamples, function(i) {
        create_spliced_count_table(og_bams[i], rot_bams[i])
    })

    # Create a consensus sequence for the isolate

    consensus <- create_consensus(piles)

    # Number of bases in the consensus
    bp <- length(consensus)

    # Calculate mutant allele frequencies at each locus for all samples
    mut_consensus <- create_mut_consensus(piles, consensus)
    mut_afs <- lapply(piles, get_mut_af, mut_consensus)
    mut_afs <- matrix(unlist(mut_afs), ncol=nsamples)

    max_sample <- get_max_sample(mut_afs)

    sb <- calc_strand_bias(og_bams, rot_bams, consensus, mut_consensus,
                           max_sample)

    # Create a two column vector for each sample with mutant allele counts vs
    # coverage at each position.
    mut_wt_counts <- lapply(piles, convert_to_mut_wt_counts, mut_consensus)
    # Convert to separate matrix for each position
    mut_wt_matrices <- lapply(1:bp, function(i) {
        mat <- sapply(mut_wt_counts, function(x) {
            x[i,]
        })
        t(mat)
    })

    # Coverage at each genome position for each sample
    coverage <- lapply(piles, get_coverage)

    # Mean coverage across samples at each locus
    mean_coverage <- sapply(c(1:bp), function(i) {
        mean(sapply(coverage, "[[", i))
    })

    # Potential heteroplasmic loci
    unique <- determine_unique(piles, mut_consensus, seq_err)

    # Low coverage variants
    overall_coverage <- mean(mean_coverage)
    coverage_proportion <- mean_coverage/overall_coverage

    samples <- sapply(og_bams, get_sample)
    # Get sample IDs for each variant
    var_samples <- sapply(c(1:bp), function(i) {
        samples[max_sample[i]]
    })

    # Get variant allele frequencies
    var_afs <- sapply(1:bp, function(i) {
        mut_afs[i, max_sample[i]]
    })

    # Get difference between mutant sample and non-mutant sample allele
    # frequencies
    diff_afs <- sapply(1:bp, function(i) {
        var_afs[i] - mean(mut_afs[i,-max_sample[i]])
    })

    # Get class of mutation
    mut_class <- sapply(mut_consensus, function(allele) {
        if (allele == "-") {
            "deletion"
        } else if (allele == "+") {
            "insertion"
        } else {
            "snv"
        }
    })

    species <- unique(get_species(og_bams))
    isolate <- unique(get_isolate(og_bams))
    genotype <- unique(sapply(og_bams, get_genotype))
    # Return output table
    test_table <- cbind(as.data.frame(rep(species, bp), stringsAsFactors=FALSE),
                        rep(genotype, bp), isolate, var_samples, c(1:bp),
                        mean_coverage, consensus, mut_consensus,
                        mut_class, var_afs, diff_afs,
                        sb, unique, coverage_proportion,
                        stringsAsFactors=FALSE)
    colnames(test_table) <- c("species", "genotype", "isolate", "sample", "pos",
                              "coverage", "ref", "alt", "class", "af", "af_diff",
                              "strand_bias", "unique", "coverage_proportion")
    list(piles, mut_wt_matrices, test_table)
}

# Merge rows of test_table given by indel_idx into a single row.
# Given:
#   test_table - A table generated by construct_test_tables
#   indel_idx - A numeric vector (length > 1)
# Return:
#   A test_table with indel_idx rows merged into single row.
#' @export
merge_indel <- function(test_table, indel_idx) {
    pre_row <- test_table[indel_idx[1] - 1, ]
    indel_row <- test_table[indel_idx,]

    # Determine alternative and reference alleles based upon indel class
    in_or_del <- unique(indel_row$alt)
    pre_ref <- pre_row$ref
    if (in_or_del == "-") {
        class <- "deletion"
        del_ref <- paste0(indel_row$ref, collapse = "")
        ref <- paste0(pre_ref, del_ref, collapse = "")
        alt <- pre_ref
    } else if (in_or_del == "+") {
        class <- "insertion"
        alt <- paste0(pre_ref, "+", collapse = "")
    } else {
        stop("These indices don't correspond to an indel.")
    }

    # The new merged row.
    new_row <- list(
        # species
        unique(as.character(indel_row$species)),
        # genotype
        unique(as.character(indel_row$genotype)),
        # isolate
        unique(as.character(indel_row$isolate)),
        # mutant sample
        unique(as.character(indel_row$sample))[1],
        # position
        pre_row$pos,
        # coverage
        mean(indel_row$coverage),
        # reference allele
        ref,
        # variant allele
        alt,
        # mutation class
        class,
        # mutant allele frequency
        mean(indel_row$af),
        # mutant sample allele frequency - mean allele frequency
        mean(indel_row$af_diff),
        # strand bias
        mean(indel_row$strand_bias),
        # is variant unique?
        all(indel_row$unique),
        # is variant low coverage?
        mean(indel_row$coverage_proportion),
        # p value
        mean(indel_row$p_value))

    new_test_table <- replace_rows(test_table, indel_idx, new_row)
    new_test_table
}

#' @export
partition_by_runs <- function(x) {
    split(x, cumsum(seq_along(x) %in%
                        (which(diff(x) > 1) + 1)))
}

# Find runs of identical values
# Given:
#   idx - A numeric vector
# Return:
#   A list of numeric vectors of values in idx. Each list item is a run.
#' @export
find_runs <- function(x) {
    # Partition idx by runs
    par <- partition_by_runs(x)
    # Return runs > 1
    runs <- par[lengths(par) > 1]
    runs
}

# Find a multinucleotide indel in a test table
# Given:
#   test_table - A table generated by construct_test_tables
# Return:
#   The indices of an arbitrary multinucleotide indel in test_table. If no
#   indels found, NULL is returned.
#' @export
find_indel <- function(test_table) {

    # First look for insertions
    is_ins <- test_table[, "alt"] == "+"
    is_sig <- test_table$p_value < 0.05

    insertion_idx <- which(is_ins & is_sig)
    insertion_run <- find_runs(insertion_idx)

    if (length(insertion_run) == 0) {
        is_del <- test_table[, "alt"] == "+"

        # If we've found all insertions, look for deletions
        deletion_idx <- which(is_del & is_sig)
        deletion_run <- find_runs(deletion_idx)

        # If all indels have been merged, return NULL
        if (length(deletion_run) == 0) {
            NULL
        } else {
            deletion_run[[1]]
        }
    } else {
        insertion_run[[1]]
    }
}

## Find runs of significant (p < 0.05) insertions or deletions and merge them
## into single indel events. Deletions are renamed to correct VCF style format,
## but insertions will need to be renamed manually.
## Given:
##   test_table - A table generated by construct_test_tables
## Return:
##   A test_table
#' @export
merge_significant_indels <- function(test_table) {
    # Find an indel
    indel_idx <- find_indel(test_table)

    if (!is.null(indel_idx)) {
        # Recursively call merge_significant_indels, until all indels merged
        updated_table <- merge_indel(test_table, indel_idx)
        merge_significant_indels(updated_table)
    } else {
        test_table
    }
}

## RSamtools pileup labels indels with "-". This function updates deletion.
## fields in test_table with VCF formatting. This allows easy export to .vcf.
## This function assumes all indels are small. Call merge_significant_indels
## first to handle multinucleotide indels. Unfortunately Rsamtools doesn't
## specify inserted bases so these have to be renamed manually.
## Given:
##   test_table - A table generated by construct_test_tables
## Return:
##   A test_table
#' @export
rename_small_deletions <- function(test_table) {
    # Get indices of all deletions in table
    del_idx <- which(test_table$class == "deletion")
    ndel <- length(del_idx)
    # Output table
    new_t <- test_table

    for (i in 1:ndel) {
        j <- del_idx[i]
        if (test_table[j, "pos"] == 0) {
            prev_idx <- nrow[test_table]
        } else {
            prev_idx <- j - 1
        }
        new_t[j, "pos"] <- test_table[prev_idx, "pos"]
        new_t[j, "ref"] <- paste0(test_table[prev_idx, "ref"],
                                  test_table[j, "ref"])
        new_t[j, "alt"] <- test_table[j, "ref"]
    }

    new_t
}

#' @export
calc_strand_bias <- function(og_bams, rot_bams, consensus, mut_consensus,
                             max_sample) {
    nsamples <- length(og_bams)
    bp <- length(consensus)

    piles <- lapply(1:nsamples, function(i) {
        create_spliced_count_table(og_bams[i], rot_bams[i],
                                   distinguish_strands = TRUE)
    })
    cnames <- colnames(piles[[1]])
    p <- sapply(1:bp, function(i) {
        mut_allele <- mut_consensus[i]
        ref_allele <- consensus[i]
        mut_idx <- which(grepl(paste0(mut_allele, "_"), cnames))
        ref_idx <- which(grepl(paste0(ref_allele, "_"), cnames))
        max_row <- unlist(piles[[max_sample[i]]][i,])
        fisher.test(matrix(c(max_row[mut_idx], max_row[ref_idx]),
                           byrow = TRUE, ncol = 2), workspace = 2e8)$p.value

    })
    phred <- p_to_q(p)
    phred
}
