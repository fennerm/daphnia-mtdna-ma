#' Annotate a variant table with snpEff and ts/tv info
#' @param var_table A variant table
#' @param snpeff_config Path to the snpeff config file
#' @export
annotate_variant_table <- function(var_table, snpeff_config) {
  annot_var_table <- add_snp_eff_annotations(var_table, snpeff_config)
  annot_var_table <- add_ts_tv_info(annot_var_table)
  annot_var_table
}

#' Add snpEff annotations to the variant table
#' @param vcf_filename VCF file input
#' @param species Species corresponding to vcf_file
#' @param config Location of snpEff config file
#' @importFrom fen.R.util read_vcf ulapply
#' @importFrom dplyr mutate
#' @return Path to the annotated snpEFF VCF file
add_snp_eff_annotations <- function(var_table, snpeff_config) {
  species <- unique(var_table$species)

  # Convert var_table to .vcf
  vcf_filename <- variant_table_to_vcf(var_table)

  # Annotate the .vcf
  snpeff_vcf <- gsub("vcf", "annot.vcf", vcf_filename)
  system(paste("snpEff", "-c", snpeff_config, paste0("d.", species),
               vcf_filename, ">", snpeff_vcf))

  # Read the .vcf data back into R
  vcfdat <- read_vcf(snpeff_vcf)
  file.remove(snpeff_vcf)

  # Parse the annotations
  annotation_field <- vcfdat[, "INFO"]
  annotations <- strsplit(annotation_field, split = "|", fixed = TRUE)

  # Add the annotations to the variant table
  annot_var_table <- mutate(var_table,
                            effect = ulapply(annotations, "[", 2),
                            severity = ulapply(annotations, "[", 3),
                            gene = ulapply(annotations, "[", 4),
                            feature = ulapply(annotations, "[", 6),
                            coding = ulapply(annotations, "[", 8))
  annot_var_table
}

#' Convert a variant table into VCF files split by sample
#'
#' VCF files are output in a folder named 'VCF/'. Headers are in 'inst/headers/'
#' @param var_table A variant table
#' @importFrom fen.R.util write_table_with_header
variant_table_to_vcf <- function(var_table) {
  species <- unique(var_table$species)

  vcf_table <- convert_to_vcf_format(var_table, species)

  header_path <- file.path("headers", paste0(species, "_header.txt"))

  # Get VCF header
  header_file <- system.file(header_path, package = "megadaph.mtdna")
  header_text <- readChar(header_file, file.info(header_file)$size)
  header_text <- trimws(header_text)

  # Write output
  vcf_filename <- file.path("/tmp", paste0(species, ".vcf"))
  write_table_with_header(vcf_table, vcf_filename, header_text, sep = "\t",
                          quote = FALSE, row.names = FALSE, col.names = FALSE)
  vcf_filename
}


#' Convert a variant table to VCF format
#' @param var_table Variant table
#' @param species The Daphnia species ('magna', 'pulex')
#' @return A data.frame in VCF format
#' @importFrom fen.R.util p_to_q
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

#' Add transition/transversion annotation to variant table
#' @param var_table A variant table
#' @return A variant table with ts_tv column appended
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

#' Construct the info field of a VCF from a row of the variant table
#' @param var_table_row A row from a variant table
#' @return Character; A .vcf formatted info field
construct_info_field <- function(var_table_row) {
  coverage <- round(as.numeric(var_table_row["coverage"]))
  allele_frequency <- var_table_row["af"]
  info <- paste0("DP=", coverage, ";AF=", allele_frequency)
  info
}
