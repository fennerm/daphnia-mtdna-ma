#!/usr/bin/env rscript
"Annotate variants with predicted effects and ts/tv

Usage:
annotate_variants.r --config=CONF VAR_TABLE

Options:
  -c --config=CONF  snpEff config file
" -> doc

main <- function(var_table, snpeff_config) {
  megadaph.mtdna::annotate_variant_table(var_table, snpeff_config)
}

if (!interactive()) {
  opts <- docopt::docopt(doc)
  main(unlist(opts["VAR_TABLE"]), unlist(opts["config"]))
}
