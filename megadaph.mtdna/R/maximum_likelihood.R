
## Set of functions for estimating mitochondrial mutation rate and effective
## population size. This method is taken from:
##
## Haag-Liautard, C., Coffey, N., Houle, D., Lynch, M., Charlesworth, B., &
## Keightley, P. D. (2008). Direct Estimation of the Mitochondrial DNA Mutation
## Rate in Drosophila melanogaster. PLoS Biology, 6(8).
## https://doi.org/10.1371/journal.pbio.0060204
##
## Notation is kept consistent with that of the paper whenever possible.

## Returns the Wright-Fisher probability of transition from i to j with
## population size N.

#' @export
transition_probability <- function(i, j, N) {
    dbinom(j, N, i/N)
}

## Given a population size N.
## Generate an NxN wright-fisher transition probability matrix.

#' @export
generate_tmatrix <- function(N) {
    outer(0:N, 0:N, FUN = transition_probability, N)
}

## Scale the cumulative frequency distribution vector (v'), by f0.

#' @export
scale_vprime <- function(vp, f0, N) {
    vp <- vp * (1-f0)
    vp[1] <- f0 + vp[1]
    vp
}

## Calculate vprime for a given population size (N) and generation number (gen)

#' @export
generate_vprime <- function(N, gen) {
    tmatrix <- generate_tmatrix(N)
    # Calculate the eigenvectors/eigenvalues of T
    eigens <- eigen(tmatrix)
    # Generate a matrix with a column for each eigenvector of T
    P <- as.matrix(eigens$vectors)
    # Take the inverse of P
    Pinv <- solve(P)
    # Generate a diagonal matrix with eigenvalues of T along the diagonal
    D <- eigens$values

    # To exponentiate the matrix we exponentiate the eigenvalues. Dn holds the
    # intermediate powers of D.
    Dn <- D

    # D_sum holds the sum of the exponentiated eigenvalues.
    D_sum <- D
    for (i in c(1:(gen-1))) {
        # Summing the eigenvalue powers is equiv to, but faster than, repeatedly
        # squaring the TMatrix, which becomes incredibly slow for large N.
        Dn <- Dn * D
        D_sum <- D_sum + Dn
    }
    # By definition of matrix exponentiation.
    U <- P %*% diag(D_sum) %*% Pinv

    # Pinv often contains imaginary numbers which are carried over to U.
    # U will never actually contain imaginary numbers, so we just take the real
    # part.
    Re(U[2,]/sum(U[2,]))
}

## Penalize unrealistic variance values, given set of allele frequencies taken
## from sequencing replicates.

#' @export
likelihood_penalty <- function(ds, sd) {
    sum(apply(ds, 1, function(rep_d) {
        dmean <- mean(rep_d)
        sum(unlist(lapply(rep_d, function(d) {
            log(dnorm(d, mean=dmean, sd=sd))
        })))
    }))
}


#' @export
pdf <- function(di, gen, v, afs, sd) {
    pdf <- sum(unlist(lapply(2:length(v), function(i) {
        v[i] * dnorm(di, mean=afs[i], sd=sd) /
            (pnorm((1-afs[i])/sd) - pnorm((0-afs[i])/sd))
    })))
    pdf
}

#' @export
mut_likelihood <- function(d, gen, v, N, sd) {
    afs <- seq(0, 1, by=1/N)
    lik <- unlist(lapply(d, pdf, gen=gen, v=v, afs=afs, sd=sd))
    log_lik <- sum(log(lik))
    log_lik
}

#' @export
wt_likelihood <- function(v, n){
    log(v[1]) * n
}

#' @export
log_likelihood <- function(param) {
    f0 <- param[1]
    if (f0 >= 1 || f0 < 0) {
        return(-Inf)
    }
    VE <- param[2]
    if (VE <= 0) {
        return(-Inf)
    }
    N <- param[3]
    bp <- param[4]
    gen <- param[5]
    sd <- sqrt(VE)
    v <- scale_vprime(vp=vp, f0=f0, N=N)

    n <- bp - nrow(d)
    lWT <- wt_likelihood(v=v, n=n)

    lPen <- likelihood_penalty(d, sd)

    lMUT <- mut_likelihood(d=d[, 1], gen=gen, v=v, N=N, sd=sd)
    lMUT + lWT + lPen
}

#' @export
n_fixed_mle <- function(N, d, gen, bp) {
    d <<- d
    vp <<- generate_vprime(N=N, gen=gen)

    fixed <- c("N", "bp", "gen")

    maxLik::maxNM(log_likelihood,
                  start=c(f0=0, VE=1, N=N, bp=bp, gen=gen),
                  fixed=fixed)
}

#' @export
mut_mle <- function(d, gen, bp, N_max, threads) {
    if (threads > 1) {
        cl <- parallel::makeCluster(threads)
    	parallel::clusterEvalQ(cl, library(maxLik))
    	parallel::clusterEvalQ(cl, library(truncnorm))
    	parallel::clusterExport(cl, c("generate_tmatrix",
    	                    "generate_vprime",
    	                    "likelihood_penalty",
    	                    "log_likelihood",
    	                    "mut_likelihood",
    	                    "n_fixed_mle",
    	                    "scale_vprime",
    	                    "transition_probability",
    	                    "wt_likelihood"))
        parallel::parLapply(cl, 2:N_max, n_fixed_mle, d=d, gen=gen, bp=bp)
    } else {
        lapply(2:N_max, n_fixed_mle, d=d, gen=gen, bp=bp)
    }
}
