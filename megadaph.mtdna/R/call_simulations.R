## Generate simulated mutant-wildtype contingency tables. Binomial mutation
## frequency is constant for each table. This corresponds to the null hypothesis
## of extreme_max_test.
## Given:
##  ntables - Integer; number of tables
##  nsamples - Integer; Number of samples = number of rows in table
##  prob - Numeric; Binomial probability of mutation
##  coverage - Numeric; If a single number is given, all tables rows will have
##             identical coverage. If range is given, table rows will have
##             variable coverage randomly distributed with the interval
## Return:
##  A list of nsamples x 2 matrices
#' @export
simulate_contingency_tables <- function(ntables, nsamples, prob=NULL,
                                        coverage=c(1000, 10000)) {
    replicate(ntables, {

        # Set prob if unset.
        if (is.null(prob)) {
            r <- runif(1)
            prob <- matrix(rep(c(r, 1 - r), nsamples), ncol = 2, byrow = TRUE)
        } else {
            if (length(prob)==1) {
                prob <- matrix(rep(c(prob, 1 - prob), nsamples), byrow=TRUE,
                               ncol = 2)
            } else {
                prob <- matrix(c(prob, 1 - prob), ncol = 2)
            }
        }
        matrix(unlist(lapply(1:nsamples, function(i) {
            rmultinom(1, round(runif(1, coverage[1], coverage[2])),
                prob = prob[i,])
        })), ncol = 2, byrow = TRUE)

    }, simplify = FALSE)
}

#' @export
get_mean_cov <- function(mut_wt) {
    mean(ulapply(mut_wt, rowSums))
}

sim_accuracy_tabs <- function(seq_err, true_pos, ntabs, nsamples, cov) {
    nfalse_pos <- ntabs - length(true_pos)
    ctabs <- simulate_contingency_tables(ntabs, nsamples, seq_err,
                                              c(cov, cov))
    false_tabs <- ctabs[1:nfalse_pos]
    true_tabs <- ctabs[(nfalse_pos + 1):ntabs]
    true_counts <- round(true_pos * cov)
    for (i in 1:length(true_pos)) {
        true_tabs[[i]][1,] <- c(true_counts[i], cov - true_counts[i])
    }
    list(false = false_tabs, true = true_tabs)
}

true_counts <- function(split_bool, ntrue) {
    sapply(split_bool, function(x) {
        true_pos <- length(which(x)) / (ntrue / length(split_bool))
        false_neg <- 1 - true_pos
        c(true_pos, false_neg)
    })
}

#' @export
#' @importFrom qvalue qvalue
test_corrado_accuracy <- function(false_tabs, true_tabs, true_freqs,
                                  fdr.level=0.01) {
    p_false <- ulapply(false_tabs, corrado_test)
    nfalse <- length(p_false)
    p_true <- ulapply(true_tabs, corrado_test)
    ntrue <- length(p_true)
    p <- c(p_false, p_true)
    q <- qvalue(p, fdr.level = 0.01)
    q_sig <- q$significant
    q_true <- q$significant[(nfalse+1):length(p)]
    q_false <- q$significant[1:nfalse]

    false_pos <- length(which(q_false)) / nfalse
    true_neg <- 1 - false_pos

    q_split <- split(q_true, as.factor(true_freqs))
    true_stats <- true_counts(q_split, ntrue)
    list(false_pos = false_pos, true_neg = true_neg, true_stats = true_stats)
}

predicted_false_pos <- function(coverage, nsites, err_rate, cutoff) {
    pbinom(cutoff, coverage, err_rate * 2, lower.tail = FALSE) * nsites
}

calc_cutoff <- function(coverage, nsites, err_rate) {
    cut <- 0
    pred_false <- predicted_false_pos(coverage, nsites, err_rate, cut)
    while (pred_false > 1) {
        cut <- cut + 1
        pred_false <- predicted_false_pos(coverage, nsites, err_rate, cut)
    }
    cut/coverage
}

binom_call <- function(ctab, cutoff) {
    freqs <- ctab[,1] / rowSums(ctab)
    nmut <- length(which(freqs > cutoff))
    if (nmut == 1) {
        TRUE
    } else {
        FALSE
    }
}

test_binom_accuracy <- function(coverage, false_tabs, true_tabs, true_freqs,
                                err_rate) {
    nfalse <- length(false_tabs)
    ntrue <- length(true_tabs)
    bp <- ntrue + nfalse
    cutoff <- calc_cutoff(coverage, bp, err_rate)
    false_calls <- ulapply(false_tabs, binom_call, cutoff)

    true_calls <- ulapply(true_tabs, binom_call, cutoff)
    true_split <- split(true_calls, true_freqs)
    false_pos <- length(which(false_calls)) / nfalse
    true_neg <- 1 - false_pos

    true_stats <- true_counts(true_split, ntrue)
    list(false_pos = false_pos, true_neg = true_neg, true_stats = true_stats)
}

#' @export
## Actual mean coverage: 4702
seq_err_call_simulation <- function(npos=100, bp=165198, nsamples=8, cov=4700) {
    true_pos_freq <- c(0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.02, 0.05, 0.1,
                       0.2, 0.4, 0.6, 0.8, 1)
    true_pos <- rep(true_pos_freq, npos)
    seq_err_probs <- seq(0, 0.004, by = 0.0002)

    accuracy <- sapply(seq_err_probs, function(se) {
        sim_tab <- sim_accuracy_tabs(seq_err = se, true_pos = true_pos,
                          ntabs = bp, nsamples = nsamples, cov = cov)
        corrado_acc <- test_corrado_accuracy(sim_tab$false, sim_tab$true,
                                             true_pos)
        binom_acc <- test_binom_accuracy(cov, sim_tab$false, sim_tab$true,
                                         true_pos, se)
        list(corrado = corrado_acc, binom = binom_acc)
    })

}

predictive_power <- function(seq_err, mut_wt) {
    sm1 <- unlist(lapply(mut_wt, function(x) (x[1,1]/(x[1,2]+x[1,1]))))
    sm_rest <- unlist(lapply(mut_wt, function(x) mean(x[2:nrow(x), 1])/
                          (mean(x[2:nrow(x), 2]) + mean(x[2:nrow(x),1]))))
    diff_mean <- abs(sm1-sm_rest)
    diff_mean <- diff_mean[which(diff_mean < 0.1)]
    diff_constant <- abs(sm1-seq_err[1])
    diff_constant <- diff_constant[which(diff_constant < 0.1)]
    plot(diff_mean, type="l")
    plot(diff_constant, type="l")
}
