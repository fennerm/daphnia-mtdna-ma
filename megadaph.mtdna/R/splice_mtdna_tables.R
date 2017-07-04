## Given:
## og, rot   Two allele frequency tables resulting from alignment to OG and ROT
##           genomes.
## spp       Either "pulex" or "magna"
## Return:
## An allele frequency table which combines the middle sections of the OG and
## ROT files.

#' @export
splice_mtdna_data <- function(og, rot, spp) {

    # The species-specific indices for splicing the tables
    if (spp == "pulex") {
        idx <- list(c(7668:11499), c(3833:11499), c(3834:7667))
    } else if (spp == "magna") {
        idx <- list(c(7475:11211), c(3738:11211), c(3738:7474))
    }

    og_class <- class(og)[1]
    rot_class <- class(rot)[1]
    # Matrix
    if ((og_class == "matrix") && (rot_class == "matrix")) {
        spliced <- rbind(rot[idx[[1]], ], og[idx[[2]],], rot[idx[[3]],])
    # Data.table
    } else if ((og_class == "data.table") &&
               (rot_class == "data.table")) {
        require(data.table)
        spliced <- rbindlist(list(rot[idx[[1]],], og[idx[[2]],],
                                  rot[idx[[3]],]))
    # Vector
    } else if ((og_class %in% c("numeric", "integer")) &&
               (rot_class %in% c("numeric", "integer"))) {
        spliced <- c(rot[idx[[1]]], og[idx[[2]]], rot[idx[[3]]])
    }
    spliced
}
