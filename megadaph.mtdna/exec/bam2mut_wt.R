#!/usr/bin/env Rscript

min_maj_indices <- function(pile) {
    min_maj_ix <- apply(pile, 1, function(x) {
        sort.int(x, decreasing = TRUE, index.return = TRUE)$ix[1:2]
    })
    min_maj_ix
}

pile_to_maj_min <- function(piles) {
    npiles <- length(piles)
    ds <- lapply(piles, destrand)
    min_maj <- lapply(1:length(ds), function(i) {
        min_maj_ix <- min_maj_indices(ds[[i]])

        counts <- lapply(1:2, function(j) {
            ix <- min_maj_ix[j,]
            bp <- length(ix)
            stranded_ix <- stranded_indices(ix)
            stranded_counts <- ulapply(1:bp, function(k) {
                piles[[i]][k, stranded_ix[k,]]
            })
            stranded_counts <- matrix(stranded_counts, ncol = 2, byrow = TRUE)
            stranded_counts
        })
        min_maj <- do.call(cbind, counts)
        min_maj
    })
    bp <- nrow(min_maj[[1]])
    min_maj_mats <- lapply(1:bp, function(i) {
        ith_row <- ulapply(min_maj, function(m) m[i,])
        matrix(ith_row, ncol = 4, byrow=TRUE)
    })

    cat("done with this pile \n")
    min_maj_mats
}

if (!interactive()) {
    args = commandArgs(trailingOnly = TRUE)
    file_list <- list_bams(args[1])

    # Create pileups using multithreading
    max_cores <- args[2]
    n <- length(file_list)
    if (n > max_cores) {
        ncore <- max_cores
    } else {
        ncore <- n
    }
    piles <- parallel::mclapply(file_list, create_pileup,
                                distinguish_strands = TRUE,
                                mc.cores = ncore)
    saveRDS(piles, "nuc_piles.Rds")
    isolates <- get_isolate(file_list)

    split_by_isolate <- split(piles, isolates)
    n <- length(split_by_isolate)
    if (n > max_cores) {
        ncore <- max_cores
    } else {
        ncore <- n
    }
    cl <- parallel::makeCluster(ncore, outfile = "cluster.log")
    parallel::clusterEvalQ(cl, library(megadaph.mtdna))
    parallel::clusterExport(cl, c("min_maj_indices",
                                  "pile_to_maj_min"))
    mut_wt_mat <- parallel::parLapply(cl, split_by_isolate, pile_to_maj_min)
    parallel::stopCluster(cl)
    saveRDS(mut_wt_mat, "nuc_mut_wt_mats.Rds")
}
