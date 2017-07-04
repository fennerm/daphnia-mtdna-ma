## Set of functions for parsing file paths

## List all bam files in given directory. If pattern given, only return bamfiles
## containing the pattern. Returns a character vector.

#' @export
list_bams <- function(dir, pattern=NULL) {
    # List all files and directories in dir
    file_list <- list.files(dir, full.names=TRUE)
    # Exclude directories
    file_list <- file_list[which(!dir.exists(file_list))]
    # Exclude non-BAM files
    file_list <- file_list[grepl('.bam$', file_list)]
    # Extract BAMs with pattern
    if (!is.null(pattern)) {
        file_list <- file_list[grepl(pattern, file_list)]
    }
    # Return normalized paths
    sapply(file_list, normalizePath)
}

# Get sample name from filename or vector of paths.

#' @export
get_sample <- function(filepath) {
    filepath <- unlist(filepath)
    sapply(filepath, function(filepath) {
        unlist(strsplit(basename(filepath), ".", fixed=TRUE))[1]
    })
}

# Get isolate ID from filename or vector of paths.

#' @export
get_isolate <- function(filepath) {
    filepath <- unlist(filepath)
    sub("^([[:alpha:]]*).*", "\\1", get_sample(filepath))
}

#' @export
get_genotype <- function(filepath) {
    filepath <- unlist(basename(filepath))
    genotype <- substr(filepath, 1, 1)
    if (genotype == "T") {
        genotype <- "TCO"
    }
    genotype
}

# Get Daphnia species from filename or vector of paths.

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
