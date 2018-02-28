#' Get sample name from filename or vector of paths.
#' @export
get_sample <- function(filepath) {
    filepath <- unlist(filepath)
    sapply(filepath, function(filepath) {
        unlist(strsplit(basename(filepath), ".", fixed=TRUE))[1] })
}

#' Get population ID from filename or vector of paths.
#' @export
get_population <- function(filepath) {
    filepath <- unlist(basename(filepath))
    population <- substr(filepath, 1, 1)
    population <- ifelse(population == "T", "TCO", population)
    population
}

#' Get genotype ID from filename or vector of paths.
#' @export
get_genotype <- function(filepath) {
  genotype <- sub("^([[:alpha:]]*).*", "\\1", get_sample(filepath))
  # Make sure that EC and SC files are not assigned the wrong genotype
  genotype <- ifelse(get_species(filepath) == "magna", substr(genotype, 1, 2),
                     genotype)
  genotype
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
