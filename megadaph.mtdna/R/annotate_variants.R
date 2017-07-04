
#' @export
write.table_with_header <- function(x, file, header, ...){
    cat(header, '\n',  file = file)
    write.table(x, file, append = T, ...)
}

#' @export
convert_to_vcf_format <- function(var_table, species) {
    if (species=="pulex") {
        contig <- "NC_000844"
    } else if (species=="magna") {
        contig <- "NC_026914"
    }

    info <- apply(var_table, 1, function(var) {
        paste0("DP=", round(as.numeric(var["coverage"])), ";AF=", var["af"])
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

#' @export
var_table_to_vcfs <- function(var_table) {
    isolate_vars <- split(var_table, var_table$isolate)
    nisolates <- length(isolate_vars)
    species <- sapply(isolate_vars, function(iso) {
        unique(iso$species)
    })
    isolates <- sapply(isolate_vars, function(iso) {
        unique(iso$isolate)
    })
    vcf_tables <- mapply(convert_to_vcf_format, isolate_vars, species)
    sapply(1:nisolates, function(i) {
        header_file <- paste0("headers/", species[i], "_header.txt")
        header_text <- readChar(header_file, file.info(header_file)$size)
        vcf_file <- paste0("VCF/", isolates[i], ".vcf")
        write.table_with_header(vcf_tables[[i]], vcf_file, header_text,
                                sep="\t", quote=FALSE, row.names=FALSE,
                                col.names=FALSE)
        vcf_txt <- readLines(vcf_file)
        vcf_txt <- vcf_txt[sapply(vcf_txt, nchar) > 1]
        file_conn <- file(vcf_file)
        writeLines(vcf_txt, file_conn)
        close(file_conn)
        normalizePath(vcf_file)
    })

}

#' @export
run_snp_eff <- function(vcf_file, species) {
    out_file <- gsub("vcf", "annot.vcf", vcf_file)
    system(paste0("java -Xmx4g -jar ~/snpEff/snpEff.jar -c \\
                  ~/snpEff/snpEff.config d.",
                  species, " ", vcf_file, " > ", out_file))
    out_file
}

#' @export
extract_snp_eff_annotations <- function(snp_eff_vcfs, var_table) {
    merged_vcf_data <- lapply(snp_eff_vcfs, read_vcf)
    merged_vcf_data <- lapply(merged_vcf_data, function(x) {
        matrix(unlist(x), ncol=8)
    })
    merged_vcf_data <- do.call(rbind,merged_vcf_data)
    annotation_field <- merged_vcf_data[,8]
    annots <- lapply(annotation_field, function(x) {
        unlist(strsplit(x, split="|", fixed=TRUE))
    })
    effect <- sapply(annots, "[", 2)
    severity <- sapply(annots, "[", 3)
    gene <- sapply(annots, "[", 4)
    feature <- sapply(annots, "[", 6)
    coding <- sapply(annots, "[", 8)
    annot_var_table <- cbind(var_table, effect, severity, gene, feature, coding)
    annot_var_table
}

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

#' @export
annotate_variants <- function(var_table) {
    vcf_files <- unlist(var_table_to_vcfs(var_table))
    out_files <- sapply(vcf_files, function(x) {
        gsub("vcf", "annot.vcf", x)

    })
    species <- get_species(vcf_files)
    snp_eff_vcfs <- mapply(run_snp_eff, vcf_files, species)
    annot_var_table <- extract_snp_eff_annotations(snp_eff_vcfs, var_table)
    annot_var_table <- add_ts_tv_info(annot_var_table)
    annot_var_table
}
