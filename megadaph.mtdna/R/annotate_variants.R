#' Write a table to file with added header using write.table
#' @param x The table
#' @param file File path to write to
#' @param header The header
#' @export
write.table_with_header <- function(x, file, header, ...){
    cat(header, '\n',  file = file)
    write.table(x, file, append = T, ...)
}

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

#' R wrapper around snpEff
#'
#' Assumes snpEff is stored in ~/src. If
#' @param vcf_file VCF file input
#' @param species Species corresponding to vcf_file
#' @param jar Location of snpEff jar file
#' @param config Location of snpEff config file
#' @return Path to the annotated snpEFF VCF file
#' @export
run_snp_eff <- function(vcf_file, species,
                        jar="~/src/snpEff/snpEff.jar",
                        config="inst/config/snpEff.config") {
    out_file <- gsub("vcf", "annot.vcf", vcf_file)
    system(paste0("java -Xmx4g -jar ", jar, " -c ", config, " d.", species, " ",
                  vcf_file, " > ", out_file))
    out_file
}

#' Add snpEff annotations to variant table
#' @param snp_eff_vcfs list of paths to snpEff annotated VCF files
#' @param var_table Variant table
#' @return A variant table with snpEff annotations added
#' @export
add_snp_eff_annotations <- function(snp_eff_vcfs, var_table) {

    # Read each VCF and merge them into a single matrix
    merged_vcf_data <- lapply(snp_eff_vcfs, read_vcf)
    merged_vcf_data <- lapply(merged_vcf_data, function(x) {
        matrix(unlist(x), ncol = 8)
    })
    merged_vcf_data <- do.call(rbind,merged_vcf_data)

    # Parse the annotations
    annotation_field <- merged_vcf_data[,8]
    annots <- lapply(annotation_field, function(x) {
        unlist(strsplit(x, split = "|", fixed=TRUE))
    })
    effect <- sapply(annots, "[", 2)
    severity <- sapply(annots, "[", 3)
    gene <- sapply(annots, "[", 4)
    feature <- sapply(annots, "[", 6)
    coding <- sapply(annots, "[", 8)

    # Add the annotations to the variant table
    annot_var_table <- cbind(var_table, effect, severity, gene, feature, coding)

    annot_var_table
}

#' Add transition/transversion annotation to variant table
#' @param var_table A variant table
#' @return A variant table with ts_tv column appended
#' @export
add_ts_tv_info <- function(var_table) {
    ts_tv <- apply(var_table, 1, function(row) {
        if (row["class"] == "snv") {
            if (all(c(row["ref"], row["alt"]) %in% c("A", "G")) ||
                all(c(row["ref"], row["alt"]) %in% c("C", "T"))) {
                "transition"
            } else {
                "transversion"
            }
        } else {
            "NA"
        }
    })
    ts_tv_var_table <- cbind(var_table, ts_tv)
    ts_tv_var_table
}

#' Wrapper script which performs the whole annotation pipeline
#' @param var_table A variant table
#' @return An annotated variant table
#' @export
annotate_variants <- function(var_table) {
    vcf_files <- unlist(var_table_to_vcfs(var_table))
    out_files <- sapply(vcf_files, function(x) {
        gsub("vcf", "annot.vcf", x)

    })
    species <- get_species(vcf_files)
    snp_eff_vcfs <- mapply(run_snp_eff, vcf_files, species)
    annot_var_table <- add_snp_eff_annotations(snp_eff_vcfs, var_table)
    annot_var_table <- add_ts_tv_info(annot_var_table)
    annot_var_table
}
