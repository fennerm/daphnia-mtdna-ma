# Given a path to a vcf file.
# Attempt to read the vcf file as a data frame. If empty, return an empty frame.

#' @export
read_vcf <- function(path) {
    tryCatch({
        table <- read.table(path, colClasses=c("character", "numeric",
                                               "character", "character",
                                               "character", "character",
                                               "character"))
        colnames(table) <- c("CHROM", "POS","ID", "REF", "ALT", "QUAL",
                                  "FILTER","INFO")
    # If VCF contains no variants read.table throws an error. Catch it and
    # produce an empty table instead
    }, error = function(e) {
        table <- data.frame(
            CHROM=factor(),
            POS=integer(),
            ID=factor(),
            REF=factor(),
            ALT=factor(),
            QUAL=integer(),
            FILTER=factor(),
            INFO=factor())
    })
    table
}

## Given:
##  vcf     A vcf table, as returned by read_vcf
## Return:
##  A vector indices of the indel entries in the vcf table

#' @export
is_indel <- function(vcf) {
    sort(c(which(nchar(vcf$REF) > 1), which(nchar(vcf$ALT) >1)))
}
## Given:
##  vcf     A vcf table, as returned by read_vcf
## Return:
##  A vector indices of the single nucleotide variant entries in the vcf table

#' @export
is_snv <- function(vcf) {
    which((nchar(vcf$REF) == 1) & (nchar(vcf$ALT) == 1))
}

## Given:
##  vcf     A vcf table, as returned by read_vcf
##  class   Either "all", "snv" or "indel"
## Return:
##  A vector of allele frequencies from the vcf table. If class is "snv" or
##  "indel", only mutations of this type are included.

#' @export
get_vcf_afs <- function(vcf, class="all") {

    # Check if class is valid input
    class <- tolower(class)

    if (!class %in% c("all", "snv", "indel")) {
        stop("Unrecognized class.")

    } else {
        info <- vcf$INFO
        # Extract the AF information from the VCF table
        af_field <- sapply(info, function(x) unlist(strsplit(x, split=";"))[2])

        # If class is defined, filter the AF vector by class
        if (class == "snv") {
            af_field <- af_field[is_snv(vcf)]
        } else if (class == "indel") {
            af_field <- af_field[is_indel(vcf)]
        }
        # Convert AF field to numeric
        af <- as.numeric(unlist(regmatches(
            af_field,
            gregexpr("[[:digit:]]+\\.*[[:digit:]]*", af_field))))
    }
    af
}
