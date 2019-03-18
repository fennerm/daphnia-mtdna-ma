#!/usr/bin/env Rscript
"Compute final summary statistics for a mutation accumulation experiment.
This includes: mutation rate calculations, figures etc. Summaries are
output as .tsv and image files in the given output directory.

Usage:
  summarize_mutations.R --variants=TSV --genome=STRING --info=TSV --outdir=OUTDIR

Options:
  -v --variants=TSV       .csv file with variant information for multiple
                          samples
  -g --genome=STRING      One of mt/nuc
  -i --info=TSV           .csv file with MA line metadata
  -o --outdir=DIR  Output directory" -> doc

library(docopt)
library(bootr)
library(futile.logger)
library(megadaph.mtdna)
library(tidyverse)

REPS <- 10000


#' Get the unique values of a (possibly-nested) column in a tibble
get_levels <- function(tbl, by) {
  tibble(level = unique(unlist(tbl[, by])))
}


#' Save a table to file
save_table <- function(table, outdir, filename) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  output_filepath <- file.path(outdir, filename)
  write_tsv(table, output_filepath)
}

#' Save a ggplot2 object to file
save_plot <- function(plot, outdir, filename) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  output_filepath <- file.path(outdir, filename)
  ggsave(
    output_filepath,
    plot = plot,
    device = "jpg",
    width = 12,
    height = 12
  )
}

#' Assign a grouping level to its higher grouping level
#'
#' E.g assign_within_grouping("population") == "species"
#'     assign_within_grouping("sample") == "genotype"
assign_within_grouping <- function(grouping, genome) {
  if (genome == "mt") {
    hierarchy <- c("species", "population", "genotype", "sample")
  } else {
    hierarchy <- c("population", "genotype", "sample")
  }
  grouping_index <- match(grouping, hierarchy)
  if (grouping_index == 1) {
    within_grouping <- NULL
  } else {
    within_grouping <- hierarchy[grouping_index - 1]
  }
  within_grouping
}

#' Assign a species to a vector of genotypes, populations or samples
assign_species <- function(grouping) {
  map_chr(unlist(grouping), function(x) {
    if (substr(x, 1, 1) %in% c("F", "G", "I")) {
      species <- "magna"
    } else {
      species <- "pulex"
    }
    species
  })
}

save_mutation_rate_plots <- function(mutation_rates, outdir, prefix = NULL) {
  grouping <- colnames(mutation_rates)[1]
  c("log10", "unscaled") %>% map(function(yscale) {
    mu_plot <- plot_mutation_rates(mutation_rates, yscale = yscale)
    plot_filename <- paste0("mutation_rates.", yscale, ".jpg")

    if (!is.null(prefix)) {
      plot_filename <- paste0(prefix, plot_filename)
    }
    save_plot(mu_plot, outdir, plot_filename)})
}

#' Calculate mutation rates and produce boxplots
analyze_mutation_rates <- function(variant_table, outdir, within="combined") {
  grouping <- colnames(variant_table)[1]
  mutation_rates <- boot_mut_rate(variant_table, reps=REPS)
  save_table(mutation_rates, outdir, "mutation_rates.tsv")
  if (grouping == "species") {
    save_mutation_rate_plots(mutation_rates, outdir)
  } else {
    # Plot each species separately
    plot_tibble <- mutation_rates %>%
      mutate(species = assign_species(!!as.name(grouping))) %>%
      group_by(species) %>%
      nest
    for (i in 1:nrow(plot_tibble)) {
      save_mutation_rate_plots(
        plot_tibble[i, "data"][[1]][[1]],
        outdir = outdir,
        prefix = paste0(plot_tibble[i, "species"][[1]][[1]], "."))
    }
  }
  mutation_rates
}

#' Compare mutation rates by bootstrap
analyze_mutation_rate_variance <- function(variant_table, genome, outdir) {
  grouping <- colnames(variant_table)[1]
  within_grouping <- assign_within_grouping(grouping, genome)
  mutation_rate_variance <- boot_compare_mut_rates(
    variant_table,
    within = within_grouping,
    reps = REPS
  )
  filename <- "mutation_rate_variance.tsv"
  save_table(mutation_rate_variance, outdir, filename)

  mutation_rate_variance
}

#' Compare mean allele frequency by bootstrap
analyze_allele_frequencies <- function(variant_table, outdir) {
  grouping <- colnames(variant_table)[1]
  within_grouping <- assign_within_grouping(grouping, "mt")
  allele_frequency_variance <- boot_compare_mean_af(
    variant_table,
    within = within_grouping,
    reps = REPS
  )
  filename <- "allele_frequency_variance.tsv"
  save_table(allele_frequency_variance, outdir, filename)
  allele_frequency_variance
}


#' Classify mutations as transitions or transversions
classify_ts_tv <- function(ref, alt) {
  if (
    (ref %in% c("C", "G")) && (alt %in% c("C", "G")) ||
    (ref %in% c("T", "G")) && (alt %in% c("T", "G")) ||
    (ref %in% c("A", "C")) && (alt %in% c("A", "C")) ||
    (ref %in% c("A", "T")) && (alt %in% c("A", "T"))
  ) {
    ts_tv <- "transversion"
  } else if (
    (ref %in% c("A", "G")) && (alt %in% c("A", "G")) ||
    (ref %in% c("C", "T")) && (alt %in% c("C", "T"))
  ) {
    ts_tv <- "transition"
  } else {
    ts_tv <- NA
  }
  ts_tv
}


#' Add transition/transversion column to mutation table
annotate_ts_tv <- function(table_of_mutations) {
  table_of_mutations[, "ts_tv"] <- map2_chr(
    table_of_mutations$ref, table_of_mutations$alt, classify_ts_tv
  )
  table_of_mutations
}


#' Calculate num indels, num A->G etc.
analyze_mutation_spectrum <- function(variant_table, outdir) {
  grouping <- colnames(variant_table)[1]
  sub_counts <- variant_table$data %>%
    map(~extract_base_sub_spectrum(., flat = TRUE))
  sub_table <- as.tibble(do.call(rbind, sub_counts))
  sub_table <- bind_cols(variant_table[, grouping], sub_table)

  # ggplot(sub_counts, aes(x=sub_type, shape=group, color=group, y=counts)) +
  #   theme_shared() +
  #   geom_point(stat="identity", position = "dodge", size=3)

  ntransversions <- sub_counts %>% map_dbl(count_transversions)
  ntransitions <- sub_counts %>% map_dbl(count_transitions)
  sub_table <- sub_table %>%
    mutate(ntransition = ntransitions) %>%
    mutate(ntransversion = ntransversions)
  filename <- "snv_spectrum.tsv"
  save_table(sub_table, outdir, filename)
  sub_table
}


mutation_file="/home/fen-arch/daphnia-mtdna-ma.private/daphnia-mtdna-ma/pipe/output/confirmed_variants/variants.tsv"
line_info_file="/home/fen-arch/daphnia-mtdna-ma.private/daphnia-mtdna-ma/pipe/input/metadata/line_info.tsv"
genome="mt"
outdir="tmp"

main <- function(mutation_file, line_info_file, genome, outdir) {
  table_of_mutations <- read_tsv(mutation_file)
  if (!("ts_tv" %in% colnames(table_of_mutations))) {
    table_of_mutations <- annotate_ts_tv(table_of_mutations)
  }
  line_info <- read_tsv(line_info_file)
  variant_table <- VariantTable(table_of_mutations, line_info)
  if (genome == "nuc") {
    groupings <- c("genotype", "population", "sample")
  } else {
    groupings <- c("genotype", "population", "species", "sample")
  }
  for (grouping in groupings) {
    grouping_outdir <- file.path(outdir, grouping)
    grouping_data <- regroup(variant_table, grouping)

    flog.info("Analyzing data at the %s level", grouping)
    mut_rates <- analyze_mutation_rates(grouping_data, grouping_outdir)

    flog.info("Comparing mutation rates between groups...")
    if (grouping != "sample") {
      analyze_mutation_rate_variance(grouping_data, genome, grouping_outdir)
      if (genome == "mt") {
        flog.info("Comparing allele frequencies between groups...")
        analyze_allele_frequencies(grouping_data, grouping_outdir)
      }
    }

    flog.info("Calculating base sub spectrum...")
    snv_table <- filter(table_of_mutations, type == "snv")
    snv_variant_table <- VariantTable(snv_table, line_info)
    grouping_data <- regroup(snv_variant_table, grouping)
    grouping_outdir <- file.path(outdir, "spectrum")
    analyze_mutation_spectrum(grouping_data, grouping_outdir)
  }

  for (mutation_type in c("snv", "insertion", "deletion")) {
    mutation_subset <- filter(table_of_mutations, type == mutation_type)
    type_variant_table <- VariantTable(mutation_subset, line_info)
    for (grouping in groupings) {
      grouping_outdir <- file.path(outdir, "type", mutation_type)
      grouping_data <- regroup(type_variant_table, grouping)
      analyze_mutation_rates(grouping_data, grouping_outdir)
    }
  }
  flog.info("Calculating mutation rate, excluding low-frequency variants...")
  high_freq_table <- filter(table_of_mutations, af > 0.1)
  high_freq_variant_table <- VariantTable(high_freq_table, line_info)
  grouping_data <- regroup(high_freq_variant_table, "species")
  grouping_outdir <- file.path(outdir, "high_freq")
  analyze_mutation_rates(grouping_data, grouping_outdir)
}

if (!interactive()) {
  opts <- docopt(doc)
  main(
    mutation_file = unlist(opts["variants"]),
    line_info_file = unlist(opts["info"]),
    genome = unlist(opts["genome"]),
    outdir = unlist(opts["outdir"])
  )
}
