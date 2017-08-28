#' Convert a variant table to VCF format
#' @param var_table Variant table
#' @param species The Daphnia species ("magna", "pulex")
#' @return A data.frame in VCF format
#' @export
convert_to_vcf_format <- function(var_table, species) {
    if (tolower(species) == "pulex") {
        contig <- "NC_000844"
    } else if (tolower(species) == "magna") {
        contig <- "NC_026914"
    }

    # Construct info field
    info <- apply(var_table, 1, function(var) {
        depth <- round(as.numeric(var["coverage"]))
        paste0("DP=", depth, ";AF=", var["af"])
    })

    nvars <- nrow(var_table)
    vcf_table <- cbind(rep(contig, nvars),
                       var_table$pos,
                       rep(".", nvars),
                       var_table$ref,
                       var_table$alt,
                       p_to_q(var_table$p_value),
                       rep("PASS", nvars),
                       info)
    vcf_table
}

#' Convert a variant table into VCF files split by genotype
#'
#' VCF files are output in a folder named "VCF/". Headers are in "inst/headers/"
#' @param var_table A variant table
#' @importFrom fen.R.util write.table_with_header
#' @export
var_table_to_vcfs <- function(var_table) {

    # Split the table by genotype
    genotype_vars <- split(var_table, var_table$genotype)
    ngenotypes <- length(genotype_vars)

    # Get species name for each species
    species <- sapply(genotype_vars, function(geno) {
        unique(geno$species)
    })

    # Get genotypes
    genotypes <- sapply(genotype_vars, function(geno) {
        unique(geno$genotype)
    })

    # Convert each variant table to a table in VCF format
    vcf_tables <- mapply(convert_to_vcf_format, genotype_vars, species)

    # Write each VCF table to file
    sapply(1:ngenotypes, function(i) {

        # Get VCF header
        header_file <- paste0("inst/headers/", species[i], "_header.txt")
        header_text <- readChar(header_file, file.info(header_file)$size)

        # Write output
        vcf_file <- paste0("VCF/", genotypes[i], ".vcf")
        write.table_with_header(vcf_tables[[i]], vcf_file, header_text,
                                sep = "\t", quote = FALSE, row.names = FALSE,
                                col.names = FALSE)
        vcf_txt <- readLines(vcf_file)
        vcf_txt <- vcf_txt[sapply(vcf_txt, nchar) > 1]
        file_conn <- file(vcf_file)
        writeLines(vcf_txt, file_conn)
        close(file_conn)
        normalizePath(vcf_file)
    })
}
