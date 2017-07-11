#!/usr/bin/env Rscript

## Script for calculating p-values of for single extreme rows in
## mutant/wildtype contingency tables using Monte Carlo simulation.

# Permute contingency matrices with fixed row and column marginals and
# calculate the maximal allele frequency.
# Given:
#   r - Number of permutations
#   row_totals, col_totals - Row and column marginals.
# Return:
#   A vector of the maximal allele frequencies obtained in the permutations

#' @export
#' @importFrom matrixStats colMaxs
simulate_max <- function(r, row_totals, col_totals) {
    prob <- row_totals/sum(row_totals)
    # Permute column 1 counts with probability proportional to row
    # marginals.
    sim_mut_counts <- rmultinom(r, size=col_totals[1], prob=prob)
    sim_wt_counts <- rmultinom(r, size=col_totals[2], prob=prob)
    # Convert to proportions
    sim_mut_prop <- sim_mut_counts/(sim_mut_counts+sim_wt_counts)
    rm(sim_mut_counts, sim_wt_counts)
    max_prop <- colMaxs(sim_mut_prop)
    rm(sim_mut_prop)
    max_prop
}


## Custom fisher-like significance test which tests for extreme maximal values
## in a k x 2 contingency table. The maximal value is defined by the row which
## shows the largest proportion of '+' events where each row is of the form
## ('+', '-'). Column and row marginals are fixed and counts in the '+' column
## are randomly redistributed. Significance is determined by how many of the
## permutations show a maximal value as extreme as the observed.
## Given:
##  ctab - A k x 2 contingency matrix.
##  reps - Number of simulations.
## Return:
##  A p-value.

#' @export
extreme_max_test <- function(ctab, reps=10000) {

    # Calculate row and column marginals
    row_totals <- rowSums(ctab)
    col_totals <- colSums(ctab)

    if (col_totals[1] == 0) {
        1
    } else {
        # Number of samples
        nsamples <- nrow(ctab)

        # Observed mutant allele frequency
        obs <- sort(ctab[,1] / row_totals)

        # Produce simulated null distribution. We repeatedly simulate pools
        # of 10000 tables to reduce memory overhead.
        sim_max_prop <- unlist(replicate(reps/10000,
                                         simulate_max(10000, row_totals,
                                                      col_totals),
                                         simplify=FALSE))


        # Proportion of simulated data at least as extreme as observed
        sim_more_extreme <- length(which(sim_max_prop >= obs[nsamples]))

        # Calculate p-value. One is added so that p-values are never 0.
        p <- (sim_more_extreme+1) / (reps +1)
        p
    }
}

## Show that extreme_max_test generates a uniform distribution of p-values under
## the null hypothesis.
## Given:
##  reps - Number of contingency tables to simulate and test
##  nsamples - Number of rows per table
## Return:
##   Histogram of p-values. Also prints the result of a Kolmogorov–Smirnov test
##   of the p-value distribution vs. the uniform distribution.

#' @export
validate_extreme_max_test <- function(reps, nsamples) {
    # Generate null distribution
    ctabs <- simulate_contingency_tables(reps, nsamples)
    # Apply the simulation test
    ps <- sapply(ctabs, extreme_max_test)
    # Print KS test
    print(ks.test(ps, runif(reps)))
    qqplot(ps, runif(reps))
    hist(ps, breaks=50)
}


## Repeatedly simulate p-values for a set of observed mutant/wild type matrices.
## Function runs until interrupted. Small p-values (< 0.01) are simulated with
## increasing precision and saved in Rds files.
## Given:
##   mut_wt_matrices - A list of k x 2 contingency tables with mutant and wild
##                     type allele counts for k samples.
##   threads - Number of threads to use.
## Return:
##   No return. Results are saved in RDS files in current directory.

#' @export
#' @importFrom parallel makeCluster clusterEvalQ parLapply
simulate_p_values <- function(mut_wt_matrices, threads, start_reps=50000) {
    cl <- makeCluster(threads)
    clusterEvalQ(cl, library(matrixStats))
    p <- unlist(parLapply(cl,
                                    mut_wt_matrices,
                                    extreme_max_test,
                                    50000))
    max_reps <- start_reps
    while(TRUE) {
        small_p_idx <- which(p < 0.01)
        sim_p <- unlist(parLapply(cl, mut_wt_matrices[small_p_idx],
                          extreme_max_test, reps=max_reps))
        p[small_p_idx] <- sim_p
        saveRDS(p, paste0("p_", as.character(max_reps), ".Rds"))
        cat(as.character(max_reps), "\n")
        max_reps <- max_reps * 10
    }
}


#' @export
effect_of_diff_seq_err <- function(nsamples=8) {
    rnorm2 <- function(n,mean,sd) {
        r <- mean+sd*scale(rnorm(n))
        while (any(r < 0)) {
            r <- mean+sd*scale(rnorm(n))
        }
        r
    }

    sds <- seq(0, 0.001, by=1e-5)

    false_dis <- lapply(sds, function(i) {
        probs <- rnorm2(nsamples, 0.002, i)
        ctabs <- simulate_contingency_tables(1000, nsamples, probs,
                                             coverage=c(1000, 1000))
        p <- sapply(ctabs, corrado_test)
        p[which(p > 1)] <- 1
        # q <- qvalue(p, fdr.level = 0.01)$qvalues
        false_dis <- length(which(p < 0.05))
        false_dis
    })

}

## We define our own function rather than using matrixStats::logSumExp because
## we were getting underflow issues.

#' @export
log_sum_exp <- function(lx) {
    max_idx <- which.max(lx)
    logsum <- log1p(sum(exp(lx[-max_idx] - lx[max_idx]))) + lx[max_idx]
    logsum[is.nan(logsum)] <- -Inf
    logsum
}

# exact_log_sum_exp <- function(x) {
#     require(Rmpfr)
#     max_idx <- which.max(x)
#     logsum <- asNumeric(log1p(sum(mpfr(exp(x[-max_idx] - x[max_idx]), 100))) +
#       mpfr(x[max_idx], 100))
#     logsum[is.nan(logsum)] <- -Inf
#     logsum
# }


#' @export
normalize <- function(log_prob) {
    if (is.matrix(log_prob)) {
        sweep(log_prob, 1, apply(log_prob, 1, log_sum_exp))
    } else {
        log_prob - log_sum_exp(log_prob)
    }
}

#' @export
probability_proportion <- function(prob, k) {
    n <- length(prob)
    prob[k]/(sum(prob[k:n]))
}

#' @export
corrado_prob <- function(j, i, nballs, log_pik, log_1_pik) {
    index_diff <- j-i

    if (index_diff >= 0) {
        if ((nballs-j == 0) && (log_pik==0)) {
            p <- 0
        } else {
            p <- lchoose(nballs-i, index_diff) + index_diff*log_pik +
                ((nballs-j)*log_1_pik)
        }
    } else {
        p <- -Inf
    }

    p
}

#' @export
#' @importFrom matrixStats colLogSumExps
build_choose_matrix <- function(nballs) {
    dim <- nballs+1
    choose_mat <- matrix(, ncol=dim, nrow=dim)
    choose_mat[dim, ] <- c(rep(-Inf, nballs), 0)
    for (i in nballs:1) {
        shifted <- c(choose_mat[i+1, 2:dim], 0)
        row_mat <- matrix(c(shifted[1:nballs],choose_mat[i+1, 1:nballs]),
                          ncol=nballs, byrow=TRUE)
        new_row <- colLogSumExps(row_mat)
        choose_mat[i, ] <- c(new_row, 0)
    }
    choose_mat
}

#' @export
build_stochastic_matrix <- function(k, prob, choose_mat) {
    nbins <- length(prob)
    dim <- ncol(choose_mat)
    nballs <- dim-1
    pik <- probability_proportion(prob, k)
    log_pik <- log(pik)
    log_pik_inv <- log(1-pik)
    if (k == 1) {
        stoch_mat <- choose_mat[1,] + (0:nballs)*log_pik +
            (nballs:0)*log_pik_inv
        stoch_mat <- matrix(stoch_mat, ncol=dim)
        stoch_mat <- normalize(stoch_mat)
    } else if (k == nbins) {
        stoch_mat <- rep(c(rep(-Inf, nballs), 0), dim)
        stoch_mat <- matrix(stoch_mat, ncol=dim, byrow=TRUE)
    } else {
        pik_diff <- (0:nballs)*log_pik
        pik_diff_inv <- (nballs:0)*(log_pik_inv)
        pik_mat <- matrix(, ncol=dim, nrow=dim)
        for (i in 1:dim) {
            pik_mat[i,] <- pik_diff_inv + pik_diff
            pik_diff <- c(-Inf, pik_diff[1:(dim-1)])
        }
        stoch_mat <- choose_mat + pik_mat
        stoch_mat <- normalize(stoch_mat)
    }
    stoch_mat
}

#' @export
build_stochastic_matrices <- function(prob, nballs) {
    choose_mat <- build_choose_matrix(nballs)
    nbins <- length(prob)
    stoch_mats <- lapply(1:nbins, build_stochastic_matrix,
                         prob=prob, choose_mat=choose_mat)
    stoch_mats
}
# build_stochastic_matrices2 <- function(prob, nballs) {
#     nbins <- length(prob)
#     mats <- lapply(1:nbins, function(k) {
#         dim <- 0:nballs
#         col <- nballs+1
#         pik <- probability_proportion(prob, k)
#         log_pik <- log(pik)
#         log_1_pik <- log(1-pik)
#         if (k == 1) {
#             mat <- vapply(dim, function(sk) {
#                 corrado_prob(sk, 0, nballs, log_pik, log_1_pik)
#             }, numeric(1))
#             mat <- matrix(mat, ncol=col)
#             mat <- normalize(mat)
#         } else if (k == nbins) {
#             mat <- rep(c(rep(-Inf, nballs), 0), col)
#             mat <- matrix(mat, ncol=col, byrow=TRUE)
#         } else {
#             mat <- unlist(lapply(dim, function(sk) {
#                 unlist(lapply(dim, function(sk_1) {
#                     corrado_prob(sk, sk_1, nballs, log_pik, log_1_pik)
#                 }))
#             }))
#             mat <- matrix(mat, ncol=col)
#             mat <- normalize(mat)
#         }
#         mat
#     })
#     mats
# }

#' @export
#' @importFrom matrixStats colLogSumExps
logmatrix_mult <- function(mat1, mat2) {
    m <- nrow(mat1)
    n <- ncol(mat2)

    matrix(colLogSumExps(as.vector(mat1) + mat2))
}

#' @export
cull <- function(mat, cutoff) {
    mat[col(mat) - row(mat) >= cutoff] <- -Inf
    mat
}

#' @export
independent_assumption_test <- function(ctab) {
    row_totals <- rowSums(ctab)
    col_totals <- colSums(ctab)
    obs <- sort(ctab[,1] / row_totals)
    prob <- row_totals/sum(row_totals)
    obs_max <- max(obs)
    cutoff <- ceiling(obs_max * row_totals)
    p <- sum(pbinom(cutoff-1, col_totals[1], prob, lower.tail = F))
    p
}

compare_tests <- function() {
    ctabs <- simulate_contingency_tables(10000, 8, 0.002,
                                         coverage = c(100, 1000))
    sim_p <- ulapply(ctabs, extreme_max_test)
    corrado_p <- ulapply(ctabs, corrado_test)
    levin_p <- ulapply(ctabs, levin_test)
    plot(sim_p, corrado_p)
    plot(sim_p, levin_p)
}
#' @export
corrado_test <- function(ctab) {
    row_totals <- rowSums(ctab)
    col_totals <- colSums(ctab)
    nsamples <- nrow(ctab)

    # Observed mutant allele frequency
    obs <- sort(ctab[,1] / row_totals)
    prob <- row_totals/sum(row_totals)
    obs_max <- max(obs)
    cutoff <- ceiling(obs_max * row_totals)
    p_ind <- sum(pbinom(cutoff-1, col_totals[1], prob, lower.tail = F))

    if (col_totals[1] == 0) {
        p <- 1
    } else if (sum(sort(cutoff)[1:2]) > col_totals[1]) {
        if (p_ind <= 0) {
            p <- .Machine$double.xmin
        } else {
            p <- p_ind
        }
    } else if (p_ind < 2.2e-16) {
        p <- 2.2e-16
    } else {

        choose_mat <- build_choose_matrix(col_totals[1])
        cum_product <- build_stochastic_matrix(1, prob, choose_mat)
        cum_product <- cull(cum_product, cutoff[1])
        for (i in 2:nsamples) {
            M <- build_stochastic_matrix(i, prob, choose_mat)
            culled <- cull(M, cutoff[i])
            rm(M)
            if (i == nsamples) {
                culled <- matrix(culled[, ncol(culled)])
            }
            cum_product <- logmatrix_mult(cum_product, culled)
            rm(culled)
        }
        rm(choose_mat)
        p <- -expm1(cum_product)
    }
    if (p <= 0) {
        p <- 2.2e-16
    } else if (p > 1) {
        p <- 1
    }

    p
}

# dtpois <- function(x, nballs, prob, lower, upper) {
#     dtrunc(x, lambda = nballs*prob, spec="pois", a=lower, b=upper)
# }

dtpois <- function(nballs, prob, cutoff, lower, upper) {
    d <- dtrunc(lower:upper, lambda = nballs*prob, spec = "pois",
                a = lower, b = upper)
    d[(cutoff + 1):(upper+1)] <- 0
    d
}
#' @importFrom truncdist dtrunc
plevin <- function(nballs, prob, cutoff) {
    nbins <- length(prob)

    conv <- dtpois(nballs, prob[1], cutoff[1], 0, nballs)
    if (nbins >= 2) {
        for (i in 2:nbins) {
            dens <- c(dtpois(nballs, prob[i], cutoff[i], 0, nballs),
                      rep(0, length(conv)-nballs+1))
            conv <- convolve(conv, rev(dens), type = "open")
        }
    }
    conv <- conv/sum(conv)
    conv_p <- conv[nballs+1]

    # conv_p <- sum(sapply(1:nbins, function(i) dtpois(nballs, prob[i], cutoff[i],
    #                                                  0, cutoff[i])))

    1-(sqrt(2 * pi * nballs)) *
    # 1- (factorial(nballs)/((nballs^nballs)*exp(-nballs))) *
        (prod(ppois(cutoff-1, nballs*prob))) *
        conv_p
        # sum(mapply(function(p, c) {
        #     dtpois(nballs, nballs, p, 0, c)
        # }, probs, cutoff))
}

levin_test <- function(ctab) {
    row_totals <- rowSums(ctab)
    col_totals <- colSums(ctab)
    nsamples <- nrow(ctab)

    # Observed mutant allele frequency
    obs <- sort(ctab[,1] / row_totals)
    prob <- row_totals/sum(row_totals)
    obs_max <- max(obs)
    cutoff <- ceiling(obs_max * row_totals)
    nmuts <- col_totals[1]
    if (nmuts == 0) {
        1
    } else {
        plevin(nmuts, prob, cutoff)
    }
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
