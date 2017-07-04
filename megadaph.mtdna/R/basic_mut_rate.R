# Simple mutation rate calculation.
# muts: numeric vector of detected mutation frequencies
# gen: average number of generations
# nuc: number of nucleotides surveyed (15333 * number of samples)

#' @export
mut_rate <- function(muts, gen, nuc) {
    if (nuc == 0) {
        0
    } else {
        sum(muts)/(gen*nuc)
    }
}

#https://stats.stackexchange.com/questions/136661/using-bootstrap-under-h0-to-perform-a-test-for-the-difference-of-two-means-repl
# bootstrap test for the equality of means
# Efron's and Tibshirani's - intro to the bootstrap, page 223

#' @export
boot_mut_rate_diff <- function(var_table, line_info, by, within=NULL) {
    boot_diff <- function(var_table, line_info, group1, group2) {
        var_table1 <- var_table[which(var_table[, by]==group1),]
        var_table2 <- var_table[which(var_table[, by]==group2),]
        line_info1 <- line_info[which(line_info[, by]==group1),]
        line_info2 <- line_info[which(line_info[, by]==group2),]
        af1 <- lapply(line_info1$sample, function(x) {
            var_table1[which(var_table1$sample==x), "af"]
        })
        af2 <- lapply(line_info2$sample, function(x) {
            var_table2[which(var_table2$sample==x), "af"]
        })
        mu01 <- mut_rate(unlist(af1), mean(line_info1$generations),
                         sum(line_info1$bp))
        mu02 <- mut_rate(unlist(af2), mean(line_info2$generations),
                         sum(line_info2$bp))
        obs <- abs(mu01 - mu02)
        mu_overall <- mut_rate(unlist(c(af1, af2)),
                               mean(c(line_info1$generations,
                                      line_info2$generations)),
                               sum(c(line_info1$bp), line_info2$bp))
        boot1 <- boot_mut_rate(af1, line_info1$generations,
                               line_info1$bp)
        boot1 <- boot1 - mu01 + mu_overall
        boot2 <- boot_mut_rate(af2, line_info2$generations,
                               line_info2$bp)
        boot2 <- boot2 - mu02 + mu_overall

        boot_diffs <- abs(boot1-boot2)
        p <- ((length(which(boot_diffs >= obs))) + 1)/(length(boot_diffs))
        p
    }

    diff_by <- function(var_table, line_info) {
        by_levels <- unique(line_info[, by])
        ngroups <- length(by_levels)
        d <- list()
        names <- c()
        for (i in c(1:ngroups)) {
            for (j in c(i:ngroups)) {
                if (j != i) {
                    d <- c(d, list(boot_diff(var_table, line_info,
                              by_levels[i], by_levels[j])))
                    names <- c(names, paste0(by_levels[i], by_levels[j]))
                }
            }
        }
        names(d) <- names
        d

    }

    if (!is.null(within)) {
        within_levels <- unique(line_info[, within])
        diffs <- lapply(within_levels, function(w) {
            var_table_within <- var_table[which(var_table[, within]==w), ]
            line_info_within <- line_info[which(line_info[, within]==w),]
            diff_by(var_table_within, line_info_within)
        })
        diffs <- unlist(diffs, recursive=FALSE)
    } else {
        diffs <- diff_by(var_table, line_info)
    }
    diffs
}


# boot_quantiles <- function(data, quants) {
#     bt <- replicate(10000, mean(sample(data, replace=TRUE)))
#     quantile(bt, probs=quants, type=5)
# }

#' @export
boot_mut_rate <- function(af, gen, nuc, reps=1e5) {
    nsamples <- length(gen)
    smp <- 1:nsamples
    bt <- unlist(replicate(reps, {
        boot_sample <- sample(smp, replace=TRUE)
        mut_rate(unlist(af[boot_sample]), mean(gen[boot_sample]), sum(nuc))
    }, simplify=FALSE))
    bt
}

boot_mutations <- function(af, gen, nuc, reps=1e5) {
    nsamples <- length(gen)
    smp <- 1:nsamples
    bt <- unlist(replicate(reps, {
        boot_sample <- sample(unlist(af), replace=TRUE)
        mut_rate(boot_sample, mean(gen), sum(nuc))
    }, simplify=FALSE))
    bt
}

#' @export
quantile_mut_rate <- function(af, gen, nuc, reps=1e5) {
    mu <- mut_rate(unlist(af), mean(gen), sum(nuc))
    boot <- boot_mut_rate(af, gen, nuc, reps)
    quantiles <- quantile(boot, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
    c(mu, quantiles)
}
