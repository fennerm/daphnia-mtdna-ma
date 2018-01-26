#' Convert a variant table into VCF files split by sample
#'
#' VCF files are output in a folder named 'VCF/'. Headers are in 'inst/headers/'
#' @param var_table A variant table
#' @importFrom fen.R.util write_table_with_header
#' @export
var_table_to_vcf <- function(var_table) {
  species <- unique(var_table$species)
  
  genotype <- unique(var_table$genotype)

  vcf_table <- convert_to_vcf_format(var_table, species)

  header_path <- file.path("headers", paste0(species, "_header.txt"))
  
  # Get VCF header
  header_file <- system.file(header_path, package = "megadaph.mtdna")
  header_text <- readChar(header_file, file.info(header_file)$size)
  header_text <- trimws(header_text)
  
  # Write output
  vcf_filename <- file.path("/tmp", paste0(genotype, ".vcf"))
  write_table_with_header(vcf_table, vcf_filename, header_text, sep = "\t", 
                          quote = FALSE, row.names = FALSE, col.names = FALSE)
  vcf_filename
}


#' Convert a variant table to VCF format
#' @param var_table Variant table
#' @param species The Daphnia species ('magna', 'pulex')
#' @return A data.frame in VCF format
convert_to_vcf_format <- function(var_table, species) {
  if (tolower(species) == "pulex") {
    contig <- "NC_000844"
  } else if (tolower(species) == "magna") {
    contig <- "NC_026914"
  }
  info <- apply(var_table, 1, construct_info_field)
  
  nvars <- nrow(var_table)
  vcf_table <- cbind(rep(contig, nvars), var_table$pos, rep(".", nvars), 
                     var_table$ref, var_table$alt, p_to_q(var_table$p_value), 
                     rep("PASS", nvars), info)
  vcf_table
}



#' Wrapper around snpEff
#' @param vcf_filename VCF file input
#' @param species Species corresponding to vcf_file
#' @param config Location of snpEff config file
#' @return Path to the annotated snpEFF VCF file
#' @export
run_snp_eff <- function(vcf_filename, species, snpeff_config) {
  output_filename <- gsub("vcf", "annot.vcf", vcf_filename)
  system(paste("snpEff", "-c", snpeff_config, paste0("d.", species), 
               vcf_file, ">", output_filename))
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
  merged_vcf_data <- do.call(rbind, merged_vcf_data)
  
  # Parse the annotations
  annotation_field <- merged_vcf_data[, 8]
  annots <- lapply(annotation_field, function(x) {
    unlist(strsplit(x, split = "|", fixed = TRUE))
  })
  effect <- sapply(annots, "[", 2)
  severity <- sapply(annots, "[", 3)
  gene <- sapply(annots, "[", 4)
  feature <- sapply(annots, "[", 6)
  coding <- sapply(annots, "[", 8)
  
  # Add the annotations to the variant table
  annot_var_table <- cbind(var_table, effect, severity, gene, feature, 
    coding)
  
  annot_var_table
}

#' Add transition/transversion annotation to variant table
#' @param var_table A variant table
#' @return A variant table with ts_tv column appended
#' @export
add_ts_tv_info <- function(var_table) {
  ts_tv <- apply(var_table, 1, function(row) {
    if (row["class"] == "snv") {
      if (all(c(row["ref"], row["alt"]) %in% c("A", "G")) || all(c(row["ref"], 
        row["alt"]) %in% c("C", "T"))) {
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

#' Construct the info field of a VCF from a row of the variant table
#' @param var_table_row A row from a variant table
#' @return Character; A .vcf formatted info field
construct_info_field <- function(var_table_row) {
  coverage <- round(as.numeric(var_table_row["coverage"]))
  allele_frequency <- var_table_row["af"]
  info <- paste0("DP=", coverage, ";AF=", allele_frequency)
  info
}
