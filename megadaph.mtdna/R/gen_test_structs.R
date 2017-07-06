interleave_cols <- function(mat1, mat2) {
    col <- ncol(mat1)
    interleaved <- lapply(1:col, function(i) {
        cbind(mat1[,i], mat2[,i])
    })
    interleaved <- matrix(unlist(interleaved), ncol = 2*col)
}

strandify <- function(mat) {
    n <- length(mat)
    divisor <- runif(n, 1, 3)
    fwd_mat <- round(mat / divisor)
    rev_mat <- mat - fwd_mat
    stranded_mat <- interleave_cols(fwd_mat, rev_mat)
    stranded_mat
}

gen_stranded_pile_row <- function() {
    setDT()
}

#' @export
#' @importFrom data.table setDT setnames
gen_pile <- function(r=1000, stranded=FALSE, range=c(1, 1000),
                          seq_err=0.002) {
    if (stranded) {
        col <- 12
        n <- 2
    } else if (!stranded) {
        col <- 6
        n <- 1
    }

    major <- sample(range[1]:range[2], r, replace = TRUE)
    minor <- rbinom(5*r, size = major, prob = seq_err)
    pile <- matrix(c(major, minor), ncol = 6)
    pile <- t(apply(pile, 1, sample, replace = FALSE))

    if (stranded) {
        pile <- strandify(pile)
    }

    pile <- setDT(as.data.frame(pile, stringsAsFactors = FALSE))
    if (stranded) {
        setnames(pile, c("A_+", "A_-", "C_+", "C_-", "G_+", "G_-", "T_+",
                            "T_-", "-_+", "-_-", "+_+", "+_-"))

    } else if (!stranded) {
        setnames(pile, c("A", "C", "G", "T", "-", "+"))
    }
    pile
}

#' @export
gen_pile_list <- function(n=8, r=1000, stranded=FALSE, range=c(1, 1000),
                          seq_err=0.002) {
    pile_list <- replicate(n, {
        gen_pile(r = r, stranded = stranded, range = range, seq_err = seq_err)
        }, simplify = FALSE)
}
#' #' @export
#' gen_test_table <-
