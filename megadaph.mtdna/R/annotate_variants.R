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
