#!/usr/bin/env Rscript
#!/usr/bin/env rscript
"Annotate variants with predicted effects and ts/tv

Usage:
  annotate_variants.r VAR_TABLE
" -> doc


main <- function(var_table) {
  # Convert to .vcf file intermediate
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


if (!interactive()) {
  opts <- docopt::docopt(doc)
  main(unlist(opts["VAR_TABLE"]))
}
