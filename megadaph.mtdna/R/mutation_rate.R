#' @param muts numeric; mutant allele frequencies
#' @param gen numeric; Mean generation number
#' @param nuc integer; Number of nucleotides surveyed
#' @return The mutation rate
#' @export
mutation_rate <- function(muts, gen, nuc) {
    arg_lengths <- c(length(gen), length(nuc))
    if (any(arg_lengths > 1)) {
        stop("Invalid gen or nuc value")
    } else if (any(c(gen, nuc) == 0)) {
        0
    } else {
        sum(muts) / (gen*nuc)
    }
}

#' Center a distribution at a given mean value
#' @param d numeric; the distribution
#' @param target numeric; the new mean
#' @return The centered distribution
center_distribution <- function(d, target) {
    d - mean(d) + target
}

#' Bootstrap the mutation rate statistic across samples
#' @param af list of numeric; mutant allele frequencies split by sample
#' @param gen list of numeric; mean generation number split by sample
#' @param nuc list of integer; bases surveyed split by sample
#' @param reps integer; Number of bootstrap replicates
#' @return A bootstrapped sample of mutation rates
#' @export
boot_mutation_rate <- function(af, gen, nuc, reps=1e5) {
    nsamples <- length(gen)
    smp <- 1:nsamples
    bt <- unlist(replicate(reps, {
        boot_sample <- sample(smp, replace = TRUE)
        mutation_rate(unlist(af[boot_sample]), mean(gen[boot_sample]), sum(nuc))
    }, simplify = FALSE))
    bt
}


#' Bootstrap difference in means between two groups
#'
#' #' https://stats.stackexchange.com/questions/136661/using-bootstrap-under-h0-to-perform-a-test-for-the-difference-of-two-means-repl
#' Efron's and Tibshirani's - intro to the bootstrap, page 223
#' @param var_table A variant table
#' @param line_info A line info table
#' @param group1,group2 Character; Group IDs in var_table
#' @param by Column name in var_table in which groups can be found
#' @param reps integer; Number of bootstrap replicates
#' @return A p-value
two_group_boot_mutation_rate_diff <- function(var_table, line_info, by, group1,
                                         group2, reps=1e5) {
    groups <- c(group1, group2)
    names(groups) <- groups

    # Exclude rows which are not in groups
    line_info <- line_info[which(line_info[, by] %in% groups),]
    var_table <- var_table[which(var_table[, by] %in% groups),]

    # Calculate the overall mutation rate across both groups
    mean_combined_gen <- mean(line_info$generations)
    combined_bp <- sum(line_info$bp)
    mu_overall <- mutation_rate(var_table$af, mean_combined_gen, combined_bp)

    # Construct bootstrapped mutation rate samples
    boot <- lapply(groups, function(g) {
        # Extract group from var_table and line_info
        group_var_table <- var_table[which(var_table[, by] == g),]
        group_line_info <- line_info[which(line_info[, by] == g),]

        # Calculate actual mutation rate
        af <- lapply(group_line_info$sample, function(smp) {
            group_var_table[which(group_var_table$sample == smp), "af"]
        })
        mean_gen <- mean(group_line_info$generations)
        bp <- sum(group_line_info$bp)
        mu <- mutation_rate(unlist(af), mean_gen, bp)

        # Bootstrap the mutation rate under the null
        boot <- boot_mutation_rate(af, group_line_info$generations, line_info$bp,
                              reps)
        # Center the bootstrap sample at the overall mutation rate to satisfy
        # the null hypothesis
        boot <- center_distribution(boot, mu_overall)
        list(mu = mu, boot = boot)
    })

    # Observed difference in mutation rates
    obs_diff <- abs(boot[[1]]$mu - boot[[2]]$mu)

    boot_diff <- abs(boot[[1]]$boot - boot[[2]]$boot)

    # Number of boot samples more extreme than observed / number of boot samples
    p <- ((length(which(boot_diff >= obs_diff))) + 1) / (length(boot_diff))
    p
}

#' @export
boot_mutation_rate_diff_by <- function(var_table, line_info, by) {
    by_levels <- unique(line_info[, by])
    ngroups <- length(by_levels)
    d <- list()
    names <- c()
    for (i in c(1:ngroups)) {
        for (j in c(i:ngroups)) {
            if (j != i) {
                d <- c(d, list(two_group_boot_mutation_rate_diff(var_table,
                                                                 line_info,
                          by_levels[i], by_levels[j])))
                names <- c(names, paste0(by_levels[i], by_levels[j]))
            }
        }
    }
    names(d) <- names
    d

}

#' Bootstrap test for the equality of means
#'
#' https://stats.stackexchange.com/questions/136661/using-bootstrap-under-h0-to-perform-a-test-for-the-difference-of-two-means-repl
#' Efron's and Tibshirani's - intro to the bootstrap, page 223
#' @param var_table A variant table
#' @param line_info A line info table
#' @param by Group in var_table to calculate the difference between
#' @param within Group var_table to calculate the difference within
#' @return Numeric; A bootstrap sample of differences in mean
#' @export
boot_mutation_rate_diff_within <- function(var_table, line_info, by, within=NULL) {
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


boot_mutations <- function(af, gen, nuc, reps=1e5) {
    nsamples <- length(gen)
    smp <- 1:nsamples
    bt <- unlist(replicate(reps, {
        boot_sample <- sample(unlist(af), replace=TRUE)
        mutation_rate(boot_sample, mean(gen), sum(nuc))
    }, simplify=FALSE))
    bt
}

#' @export
quantile_mutation_rate <- function(af, gen, nuc, reps=1e5) {
    mu <- mutation_rate(unlist(af), mean(gen), sum(nuc))
    boot <- boot_mutation_rate(af, gen, nuc, reps)
    quantiles <- quantile(boot, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
    c(mu, quantiles)
}
