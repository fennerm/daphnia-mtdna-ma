#' Get sample name from filename or vector of paths.
#' @export
get_sample <- function(filepath) {
    filepath <- unlist(filepath)
    sapply(filepath, function(filepath) {
        unlist(strsplit(basename(filepath), ".", fixed=TRUE))[1] })
}

#' Get isolate ID from filename or vector of paths.
#' @export
get_isolate <- function(filepath) {
    filepath <- unlist(basename(filepath))
    isolate <- substr(filepath, 1, 1)
    isolate <- ifelse(isolate == "T", "TCO", isolate)
    isolate
}

#' Get genotype ID from filename or vector of paths.
#' @export
get_genotype <- function(filepath) {
    sub("^([[:alpha:]]*).*", "\\1", get_sample(filepath))
}

#' Get Daphnia species from filename or vector of paths.
#' @export
get_species <- function(filepath) {
    filepath <- unlist(filepath)
    base <- basename(filepath)
    sapply(base, function(base) {
        if (any(startsWith(base, c("L", "T")))) {
            "pulex"
        } else if (any(startsWith(base, c("F", "G", "I")))) {
            "magna"
        }
    })
}
