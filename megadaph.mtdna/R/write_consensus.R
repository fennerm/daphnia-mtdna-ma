
#' @export
write_consensus_fastas <- function(og_bams, rot_bams, isolate_names) {
    nisolates <- length(og_bams)

    for (i in 1:nisolates) {
        nsamples <- length(og_bams[[i]])
        piles <- lapply(1:nsamples, function(j) {
            create_spliced_count_table(og_bams[[i]][j], rot_bams[[i]][j])
        })
        filename <- paste0(isolate_names[i], ".fa")
        create_consensus_fasta(piles, filename)
    }
}
