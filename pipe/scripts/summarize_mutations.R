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
  -o --output_dir=OUTDIR  Output directory" <- doc
# library(megadaph.mtdna)
# import::from(fen.R.util, save_results, select_groups, precomputed_boxplot)
import::from(docopt, docopt)
import::from(tools, file_ext)
import::from(Hmisc, capitalize)

library(tidyverse)

nvim.srcdir("lib/megadaph.mtdna/R")
nvim.srcdir("~/fmacrae/code/fen.R.util/R")
variant_table <- as.tibble(read.csv("~/fmacrae/daphnia-mtdna-ma.private/daphnia-mtdna-ma/pipe/output/merge_variants/variants.csv", stringsAsFactors = FALSE))
line_info <- as.tibble(read.csv("~/fmacrae/daphnia-mtdna-ma.private/daphnia-mtdna-ma/pipe/input/metadata/line_info.csv", stringsAsFactors = FALSE))
level <- "population"
library(boot)
outdir <- "output/summarize_mutations"
REPS <- 100

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

## Main
main <- function(variant_table, line_info, outdir) {
  combined_table <- merge_tables(variant_table, line_info)
  # levels <- c(
  #   "genotype", "population", "species", "ts_tv", "severity", "gene", "effect")

  levels <- c("genotype", "population", "species")
  inner_levels <- c("combined", "class")

  for (level in levels) {
    print(level)
    level_outdir <- file.path(outdir, level)
    level_data <- combined_table %>% group_by(!!as.name(level)) %>% nest
    for (inner_level in inner_levels) {
      print(inner_level)
      if (inner_level == "combined") {
        inner_level_outdir <- level_outdir
        inner_level_data <- level_data
      } else {
        inner_level_outdir <- file.path(outdir, level, inner_level)
        inner_level_data <- merge_tables_by_mutation_type(
          variant_table, line_info, level, inner_level)
      }
      analyze_mutation_rates(inner_level_data, inner_level_outdir)
    }
    analyze_mutation_rate_variance(level_data, level_outdir)
    analyze_allele_frequencies(level_data, level_outdir)

    # analyze_mutation_types(
    #   variant_table = variant_table,
    #   line_info = line_info,
    #   by = level,
    #   outdir = level_outdir)
  }
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
      xlab = "Mutation Rate",
      ylab = capitalize(grouping),
      yscale = yscale,
      xval = mutation_rates[, grouping],
      fill = mutation_rates[, grouping],
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

if (!interactive) {
  opts <- docopt(doc)
  variant_table <- read.csv(unlist(opts["variants"]), stringsAsFactors = FALSE)
  line_info_table <- read.csv(unlist(opts["info"]), stringsAsFactors = FALSE)
  main(
    variant_table, line_info_table, unlist(opts["output_dir"]))
}


#' @export
by_effect <- function(variant_table, by) {

  split_by <- split(variant_table, variant_table[, by])

  ts_tv_by <- t(sapply(split_by, function(t) {
      nts <- length(which(t$ts_tv=="transition"))
      ntv <- length(which(t$ts_tv=="transversion"))
      c(nts, ntv)
      }))

  nts <- length(which(variant_table$ts_tv=="transition"))
  ntv <- length(which(variant_table$ts_tv=="transversions"))
  prop.test(nts, nts+ntv, p=1/3)

  a_t <- length(which((variant_table$ref %in% c("A", "T")) &
      (variant_table$alt %in% c("A", "T"))))
  g_t <- length(which((variant_table$ref %in% c("G", "T")) &
      (variant_table$alt %in% c("G", "T"))))
  g_c <- length(which((variant_table$ref %in% c("G", "C")) &
      (variant_table$alt %in% c("G", "C"))))
  a_c <- length(which((variant_table$ref %in% c("A", "C")) &
      (variant_table$alt %in% c("A", "C"))))
  a_g <- length(which((variant_table$ref %in% c("A", "G")) &
      (variant_table$alt %in% c("A", "G"))))
  c_t <- length(which((variant_table$ref %in% c("C", "T")) &
      (variant_table$alt %in% c("C", "T"))))

  snv_df <- data.frame(id=c("A->T", "G->T", "G->C", "A->C", "A->G",
      "C->T"), n=c(a_t, g_t, g_c, a_c, a_g, c_t),
    stringsAsFactors = FALSE)
  snv_df <- snv_df[order(snv_df$n),]
}

#' @export
#' @importFrom coin independence_test
by_mutation_effect <- function(variant_table) {
  independence_test(variant_table$af~as.factor(variant_table$effect))
  independence_test(variant_table$af~as.factor(variant_table$gene))
  synon <- as.factor(variant_table$effect=="synonymous_variant")
  independence_test(variant_table$af~synon)
  independence_test(variant_table$p_value~synon)
}

#' @export
locations <- function(variant_table, n) {
  mag_mut_by_loc <- sapply(1:(n-1), function(i) {
    nrow(dplyr::filter(variant_table, species=="magna", pos > (i*(14948/n)),
        pos <= ((i+1)*(14948/n))))
    })
  pul_mut_by_loc <- sapply(1:(n-1), function(i) {
    nrow(dplyr::filter(variant_table, species=="pulex", pos > (i*(15333/n)),
        pos <= ((i+1)*(15333/n))))
    })
}

#' @export
read_tables <- function() {
  variant_table <<- read.csv("../data/tables/annot_table.csv", stringsAsFactors = F)
  line_info <<- read.csv("../data/tables/line_info.csv", stringsAsFactors = F)
  merged_table <<- read.csv("../data/tables/merged_table.csv",
    stringsAsFactors = F)
}

#' @export
species_mu <- function(line_info, variant_table) {
  species_mutation_rates <- mutation_rate_by(variant_table, line_info, "species")
  box <- custom_boxplot(species_mutation_rates, "Species", "Mutation Rate",
    c("D. magna", "D. pulex"))
  box

}

#' @export
species_diffs <- function(line_info, variant_table) {
  species_diffs <- boot_mutation_rate_diff_within(variant_table, line_info, "species")
  indel_diffs <- boot_mutation_rate_diff_within(dplyr::filter(variant_table, class != "snv"),
    line_info, "species")

  c(species_diffs, indel_diffs)
}
