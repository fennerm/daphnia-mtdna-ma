npairs <- 250000
npos <- 15000
n <- npairs / npos
mean_insert <- 200
sttdev_insert <- 55

simulate_position <- function(n, mean_insert, stddev_insert) {
    sum(table(round(rnorm(n, mean_insert, stddev_insert))) - 1)
}

sum(replicate(npos, simulate_position(n, mean_insert, sttdev_insert)))
