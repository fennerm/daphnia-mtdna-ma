#!/usr/bin/env Rscript

## Given the species of Daphnia and the position on the rotated reference seq
## Return the position of this base on the original reference sequence.

#' @export
rot_to_og <- function(pos, spp) {
    if (tolower(spp) == "pulex") {
        if (pos < 7668) {
            pos + 7666
        } else if (pos > 7667 & pos < 15334) {
            pos - 7667
        }
    } else if (tolower(spp) == "magna") {
        if (pos < 7475) {
            pos + 7474
        } else if (pos > 7474 & pos < 14949) {
            pos - 7474
        }
    }
}

## Given the species of Daphnia and the position on the original reference seq
## Return the position of this base on the rotated sequence.

#' @export
og_to_rot <- function(pos, spp) {
    if (tolower(spp) == "pulex") {
        if (pos < 7667) {
            pos + 7667
        } else if (pos > 7668 & pos < 15334) {
            pos - 7666
        }
    } else if (tolower(spp) == "magna") {
        if (pos < 7474) {
            pos + 7475
        } else if (pos > 7474 & pos < 14949) {
            pos - 7474
        }
    }
}
