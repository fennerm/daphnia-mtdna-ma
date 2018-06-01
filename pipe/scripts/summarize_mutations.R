#!/usr/bin/env Rscript
"Compute final summary statistics for the Daphnia mtDNA mutation accumulation
project. This includes: mutation rate calculations, figures etc. Summaries are
output as .csv and image files in the given output directory.

Usage:
  summarize_mutations.R --variants=CSV --info=CSV --output_dir=OUTDIR

Options:
  -v --variants=CSV       .csv file with variant information for multiple
                          samples
  -i --info=CSV           .csv file with MA line metadata
  -o --output_dir=OUTDIR  Output directory" -> doc


# library(megadaph.mtdna)
# import::from(fen.R.util, save_results, select_groups, precomputed_boxplot)
import::from(docopt, docopt)
import::from(tools, file_ext)
import::from(Hmisc, capitalize)

library(tidyverse)
library(futile.logger)

library(fen.R.util)
library(megadaph.mtdna)

REPS <- 1000
variant_table <- as.tibble(read.csv("~/fmacrae/daphnia-mtdna-ma.private/daphnia-mtdna-ma/pipe/output/annotate_variants/variants.csv", stringsAsFactors = FALSE))
line_info <- as.tibble(read.csv("~/fmacrae/daphnia-mtdna-ma.private/daphnia-mtdna-ma/pipe/input/metadata/line_info.csv", stringsAsFactors = FALSE))

#' Merge the variant and line info tables into a table
merge_tables <- function(variant_table, line_info) {
  variant_table <- variant_table %>%
      select(-species, -population, -genotype) %>%
      group_by(sample) %>%
      summarize_all("list")
  line_info <- line_info %>% rename(mean_coverage = coverage)
  mutation_table <- left_join(line_info, variant_table, by = "sample")
  mutation_table$af <- mutation_table$af %>% replace_null(0)
  mutation_table
}

#' Get the unique values of a (possibly-nested) column in a tibble
get_levels <- function(tbl, by) {
  tibble(level = unique(unlist(tbl[, by])))
}

#' Merge line_info and variant_tables seperately for each level of a mutation
#' type
#'
#' E.g If level = "species" and inner_level = "class". This function returns a
#' tibble with columns 1:2 =
#' ["magna", "pulex"] x ["insertion", "deletion", "snv"] and column 3 the merged
#' table for the combination.
merge_tables_by_mutation_type <- function(
    variant_table,
    line_info,
    level,
    inner_level
  ) {
  merged_table <- get_levels(line_info, level) %>%
    # Split the variant table by `level`
    mutate(table = level %>%
      map(~filter(variant_table, !!as.name(level) == .))) %>%
    # Create a merged mutation table for each variant table separately
    mutate(data = table %>% map(function(x) {
      level_data <- as.tibble(unique(variant_table[, inner_level])) %>%
        mutate(intermediate = !!as.name(inner_level) %>%
          map(~filter(x, !!as.name(inner_level) == .))) %>%
        mutate(data = intermediate %>%
          map(~merge_tables(., line_info))) %>%
        select(class, data)
    })) %>%
    # Remove the intermediate results
    select(level, data) %>%
    unnest
  colnames(merged_table)[1] <- level
  merged_table
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
    height = 12)
}

#' Assign a grouping level to its higher grouping level
#'
#' E.g assign_within_grouping("population") == "species"
#'     assign_within_grouping("sample") == "genotype"
assign_within_grouping <- function(grouping) {
  hierarchy <- c("species", "population", "genotype", "sample")
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
    species})
}

plot_mutation_rates <- function(mutation_rates, outdir, prefix = NULL) {
  grouping <- colnames(mutation_rates)[1]
  c("log10", "unscaled") %>% map(function(yscale) {
    plot <- precomputed_boxplot(
      mutation_rates,
      ylab = "Mutation Rate",
      xlab = capitalize(grouping),
      yscale = yscale,
      xval = mutation_rates[, grouping],
      fill = mutation_rates[, grouping],
      limits = c(8e-9, 8e-6),
      breaks = c(1e-8, 1e-7, 1e-6, 5e-6, 1e-5),
      legend = FALSE
    )
    plot_filename <- paste0("mutation_rates.", yscale, ".jpg")
    if (!is.null(prefix)) {
      plot_filename <- paste0(prefix, plot_filename)
    }
    save_plot(plot, outdir, plot_filename)})
}

#' Calculate mutation rates and produce boxplots
analyze_mutation_rates <- function(mutation_table, outdir) {
  grouping <- colnames(mutation_table)[1]
  mutation_rates <- boot_quantiles(
    mutation_table, indexed_mutation_rate, reps = REPS)

  mutation_rates_filename <- "mutation_rates.tsv"
  save_table(mutation_rates, outdir, mutation_rates_filename)

  if (grouping == "species") {
    plot_mutation_rates(mutation_rates, outdir)
  } else {
    # Plot each species separately
    plot_tibble <- mutation_rates %>%
      mutate(species = assign_species(!!as.name(grouping))) %>%
      group_by(species) %>%
      nest

    for (i in 1:nrow(plot_tibble)) {
      print(i)
      plot_mutation_rates(
        plot_tibble[i, "data"][[1]][[1]],
        outdir = outdir,
        prefix = paste0(plot_tibble[i, "species"][[1]][[1]], "."))
    }
  }
  mutation_rates
}

#' Compare mutation rates by bootstrap
analyze_mutation_rate_variance <- function(mutation_table, outdir) {
  grouping <- colnames(mutation_table)[1]
  within_grouping <- assign_within_grouping(grouping)
  mutation_rate_variance <- boot_compare_all(
    mutation_table,
    within = within_grouping,
    statistic = indexed_mutation_rate,
    reps = REPS
  )

  filename <- "mutation_rate_variance.tsv"
  save_table(mutation_rate_variance, outdir, filename)

  mutation_rate_variance
}

#' Compare mean allele frequency by bootstrap
analyze_allele_frequencies <- function(mutation_table, outdir) {
  grouping <- colnames(mutation_table)[1]
  within_grouping <- assign_within_grouping(grouping)
  allele_frequency_variance <- boot_compare_all(
    mutation_table,
    within = within_grouping,
    statistic = indexed_mean_af,
    reps = REPS
  )
  filename <- "allele_frequency_variance.tsv"
  save_table(allele_frequency_variance, outdir, filename)
  allele_frequency_variance
}

#' Calculate num indels, num A->G etc.
analyze_mutation_types <- function(mutation_table, outdir) {
  grouping <- colnames(mutation_table)[1]

  sub_counts <- mutation_table$data %>%
    map(~extract_base_sub_spectrum(.)) %>%
    do.call(rbind, .)
  sub_table <- as.tibble(cbind(mutation_table[1], sub_counts))

  ts_tv <- mutation_table$data %>%
    map(~unlist(.$ts_tv))
  mut_class <- mutation_table$data %>%
    map(~unlist(.$class))

  sub_table <- sub_table %>%
    mutate(ntransition = ts_tv %>%
      map_int(count_occurences, "transition")) %>%
    mutate(ntransversion = ts_tv %>%
      map_int(count_occurences, "transversion")) %>%
    mutate(nsnv = mut_class %>%
      map_int(count_occurences, "snv")) %>%
    mutate(ninsertion = mut_class %>%
      map_int(count_occurences, "insertion")) %>%
    mutate(ndeletion = mut_class %>%
      map_int(count_occurences, "deletion"))

  save_table(sub_table, outdir, "mutation_counts.tsv")
  # nts <- length(which(variant_table$ts_tv=="transition"))
  # ntv <- length(which(variant_table$ts_tv=="transversions"))
  # prop.test(nts, nts+ntv, p=1/3)
  #
  # snv_df <- data.frame(id=c("A->T", "G->T", "G->C", "A->C", "A->G",
  #     "C->T"), n=c(a_t, g_t, g_c, a_c, a_g, c_t),
  #   stringsAsFactors = FALSE)
  # snv_df <- snv_df[order(snv_df$n),]
}


## Main
main <- function(variant_table, line_info, outdir) {
  combined_table <- merge_tables(variant_table, line_info)
  # levels <- c(
  #   "genotype", "population", "species", "ts_tv", "severity", "gene", "effect")

  levels <- c("genotype", "population", "species")
  inner_levels <- c("combined", "class")

  for (level in levels) {
    flog.info("Analyzing data at the %s level", level)
    level_outdir <- file.path(outdir, level)
    level_data <- combined_table %>% group_by(!!as.name(level)) %>% nest
    for (inner_level in inner_levels) {
      flog.info("Analyzing mutation rates by %s", inner_level)
      if (inner_level == "combined") {
        inner_level_outdir <- level_outdir
        inner_level_data <- level_data
      } else {
        inner_level_outdir <- file.path(outdir, level, inner_level)
        inner_level_data <- merge_tables_by_mutation_type(
          variant_table, line_info, level, inner_level)
      }
      # TODO - Fix by class figure output
      analyze_mutation_rates(inner_level_data, inner_level_outdir)
    }
    flog.info("Comparing mutation rates between groups...")
    analyze_mutation_rate_variance(level_data, level_outdir)
    flog.info("Comparing allele frequencies between groups...")
    analyze_allele_frequencies(level_data, level_outdir)
    # flog.info("Counting mutation types...")
    # analyze_mutation_types(level_data, level_outdir)
  }
}

if (!interactive()) {
  opts <- docopt(doc)
  main(
    variant_table = read_csv(unlist(opts["variants"])),
    line_info = read_csv(unlist(opts["info"])),
    outdir = unlist(opts["output_dir"]))
}


