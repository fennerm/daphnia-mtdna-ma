interleave.cols <- function(mat1, mat2) {
    col <- ncol(mat1)
    interleaved <- lapply(1:col, function(i) {
        cbind(mat1[,i], mat2[,i])
    })
    interleaved <- matrix(unlist(interleaved), ncol = 2*col)
}

strandify <- function(mat) {
    n <- length(mat)
    divisor <- runif(n, 1, 3)
    fwd.mat <- round(mat / divisor)
    rev.mat <- mat - fwd.mat
    stranded.mat <- interleave.cols(fwd.mat, rev.mat)
    stranded.mat
}

#' @export
#' @importFrom data.table setDT setnames
gen.stranded.pile.row <- function() {
    r <- setDT(list(200, 100, 0, 0, 10, 15, 0, 0, 100, 50, 0, 0))
    setnames(r, c("A.+", "A.-", "C.+", "C.-", "G.+", "G.-", "T.+",
                  "T.-", "-.+", "-.-", "+.+", "+.-"))
    r
}

#' @export
#' @importFrom data.table setDT setnames
gen.unstranded.pile.row <- function() {
    r <- setDT(list(300, 0, 25, 0, 150, 0))
    setnames(r, c("A", "C", "G", "T", "-", "+"))
    r
}

#' @export
#' @importFrom data.table setDT setnames
gen.pile <- function(r=100, stranded=FALSE, range=c(1, 1000),
                          seq.err=0.002) {
    if (stranded) {
        col <- 12
        n <- 2
    } else if (!stranded) {
        col <- 6
        n <- 1
    }

    major <- sample(range[1]:range[2], r, replace = TRUE)
    minor <- rbinom(5*r, size = major, prob = seq.err)
    pile <- matrix(c(major, minor), ncol = 6)
    pile <- t(apply(pile, 1, sample, replace = FALSE))

    if (stranded) {
        pile <- strandify(pile)
    }

    pile <- setDT(as.data.frame(pile, stringsAsFactors = FALSE))
    if (stranded) {
        setNames(pile, c("A.+", "A.-", "C.+", "C.-", "G.+", "G.-", "T.+",
                            "T.-", "-.+", "-.-", "+.+", "+.-"))

    } else if (!stranded) {
        setNames(pile, c("A", "C", "G", "T", "-", "+"))
    }
    pile
}

#' @export
gen.pile.list <- function(n=8, r=1000, stranded=FALSE, range=c(1, 1000),
                          seq.err=0.002) {
    pile.list <- replicate(n, {
        gen.pile(r = r, stranded = stranded, range = range, seq.err = seq.err)
        }, simplify = FALSE)
}
