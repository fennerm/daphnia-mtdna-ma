# Convert phred score to probability of error.

#' @export
q_to_p <- function(q) {
    10 ^ ((-q) / 10)
}

# Convert probability of error to phred score

#' @export
p_to_q <- function(p) {
    q <- -10 * log10(p)
    q <- sapply(q, function(x) {
        if (is.infinite(x)) {
            100000
        } else {
            x
        }
    })
    q
}
