#!/usr/bin/env Rscript

#' Should be a script
#' @param og_bams List of bam files from alignment to the unrotated reference
#'                sequences
#' @param rot_bams List of bam files from alignment to the rotated reference
#'                sequences
#' @param seq_err_rates List of sequencing error rate estimates for each sample
#' @export
call_variants <- function(og_bams, rot_bams, seq_err_rates) {
    nsamples <- length(og_bams)

    tabs <<- lapply(1:nsamples, function(i) {
        construct_tables(og_bams[[i]], rot_bams[[i]],
                         seq_err_rates[[i]], subsample)
    })
    piles <<- lapply(tabs, "[[", 1)
    mut_wt_matrices <<- lapply(tabs, "[[", 2)
    test_tables <<- lapply(tabs, "[[", 3)
    merged_table <- do.call(rbind, test_tables)
    p_value <- readRDS("../tables/corrado_p.Rds")
    merged_table <- cbind(merged_table, p_value)
    test_tables <- split(merged_table, merged_table$isolate)
    formatted_test_tables <- lapply(test_tables, merge_significant_indels)
    formatted_test_tables <<- lapply(formatted_test_tables, rename_small_deletions)
    fdr_table <<- mult_comparisons_correct(formatted_test_tables,
                                           fdr_level=0.01)
    var_table <<- get_stat_significant(fdr_table)
    var_table
}

# args = commandArgs(trailingOnly=TRUE)
# mut_wt_matrices <- readRDS(as.character(args[1]))
# threads <- as.numeric(args[2])
# sq <- 1:length(mut_wt_matrices)
# partitions <- split(sq, cut(sq, 40, labels=FALSE))
# cl <- makeCluster(threads, outfile="cluster.log")
# clusterExport(cl, c("logmatrix_mult", "independent_assumption_test",
#                     "normalize", "cull", "build_stochastic_matrix",
#                     "build_choose_matrix", "probability_proportion",
#                     "log_sum_exp", "corrado_test.R", "mut_wt_matrices",
#                     "partitions"))
# p <- unlist(parLapply(cl, partitions, function(is) {
#     unlist(lapply(mut_wt_matrices[is], corrado_test))
# }))
# saveRDS(p, "corrado_p.Rds")
# stopCluster(cl)

