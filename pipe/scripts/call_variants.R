#!/usr/bin/env Rscript
## Variants were called using a contingency matrix approach.
## For each group of samples with a common ancestor genotype:
## 1. Determine the major and minor allele at every locus in each sample.
##    - The major allele is defined as the nucleotide which is present in the
##      the most reads at a site.
## 2. Determine the consensus major allele at each site across all samples.
## 3. Define the minor allele as the allele which reached the highest allele
##    frequency in a single sample, other than the major allele.
## 4. For each site in the genome, produce a contingency table of the form
##              | Major | Minor |
##     Sample 1 |   -   |   -   |
##        ...   |   -   |   -   |
##     Sample N |   -   |   -   |
##
##    where 'Major' and 'Minor' are the major and minor allele frequency in
##    each sample.
## 5. Test the null hypothesis that the allele frequencies at a particular site
##    are equal across samples.
##    - For this purpose we adapted a statistical test based on the expected
##      distribution of the maximum of a multinomial distribution:
##      The exact distribution of the maximum, minimum and the range of
##      Multinomial/Dirichlet and Multivariate Hypergeometric frequencies,
##      Corrado, 2011
"Call variants from allele count pileup

Usage:
  call_variants.R --seqerr=SEQERR --output_dir=OUTDIR PILEUP ...

Options:
  -e --seqerr=SEQERR        .csv file with estimates of sequencing error for
                            each sample
  -o --output_dir=OUTDIR    Output directory
" -> doc

library(fen.R.util)
library(megadaph.mtdna)

MUT_COV_MATRICES_DIR <- "mut_cov_matrices"
TEST_TABLE_DIR <- "test_tables"
CONSENSUS_DIR <- "consensus_seqs"

## Write data to an .Rds or .csv file
write_output <- function(dat, dir, genotype, type) {
  if (toupper(type) == "RDS") {
        filename <- paste0(genotype, ".Rds")
        saveRDS(dat, file.path(dir, filename))
  } else if (toupper(type) == "CSV") {
    filename <- paste0(genotype, ".csv")
    write.csv(dat, file = file.path(dir, filename), row.names = FALSE,
              quote = FALSE)
  }
}

## Create the output directories
create_output_directories <- function(outdir) {
  dir.create(file.path(outdir, MUT_COV_MATRICES_DIR))
  dir.create(file.path(outdir, TEST_TABLE_DIR))
  dir.create(file.path(outdir, CONSENSUS_DIR))
}

## Write the consensus sequence to file
write_consensus <- function(consensus, outdir, genotype) {
  filename <- file.path(outdir, CONSENSUS_DIR, paste0(genotype, ".fa"))
  consensus_string <- paste0(consensus, collapse = "")
  fileconn <- file(filename)
  writeLines(c(paste0(">", genotype), consensus_string), fileconn)
  close(fileconn)
}

## Read the sequencing error file and convert to a numeric vector
parse_seq_err_csv <- function(seq_err_file) {
  seq_err <- read.csv(seq_err_file, stringsAsFactors = FALSE)
  seq_err <- seq_err$sequencing_error
  seq_err
}

main <- function(pileup_files, seq_err_file, outdir) {
  nsamples <- length(pileup_files)
  samples <- get_sample(pileup_files)
  isolate <- unique(get_isolate(pileup_files))
  genotype <- unique(get_genotype(pileup_files))
  species <- unique(get_species(pileup_files))

  cat("Creating output directories \n")
  create_output_directories(outdir)

  cat("Reading pileups \n")
  stranded_piles <- lapply(pileup_files, read.csv, stringsAsFactors = FALSE,
                           check.names = FALSE)
  stranded_piles <- lapply(stranded_piles, data.table::setDT)
  piles <- lapply(stranded_piles, destrand)

  cat("Reading sequencing error file \n")
  seq_err <- parse_seq_err_csv(seq_err_file)

  cat("Building consensus sequence \n")
  consensus <- create_consensus(piles)
  bp <- length(consensus)
  write_consensus(consensus, outdir, genotype)

  cat("Determining potential mutant alleles \n")
  mutant_consensus <- create_mutant_consensus(piles, consensus)

  cat("Calculating allele count statistics \n")
  mutant_counts <- lapply(piles, get_allele_counts, mutant_consensus)
  coverage <- lapply(piles, get_coverage)
  mean_coverage <- unlist(apply_across(coverage, mean))
  mut_cov_counts <- mapply(cbind, mutant_counts, coverage, SIMPLIFY = FALSE)
  mut_cov_matrices <- by_position(mut_cov_counts)

  write_output(mut_cov_matrices, dir = file.path(outdir, MUT_COV_MATRICES_DIR),
               genotype = genotype, type = "RDS")

  overall_coverage <- mean(mean_coverage)
  coverage_proportion <- mean_coverage / overall_coverage

  cat("Calculating mutant allele frequency \n")
  mut_afs <- mapply("/", mutant_counts, coverage)

  cat("Determining maximum mutant allele frequency per site\n")
  mutant_samples <- apply(mut_afs, 1, which.max)
  max_mut_afs <- ulapply(1:bp, function(i) mut_afs[i, mutant_samples[i]])
  diff_afs <- ulapply(1:bp, function(i) {
                          max_mut_afs[i] - mean(mut_afs[i, -mutant_samples[i]])
  })

  cat("Calculating strand bias \n")
  strand_bias <- compute_all_strand_bias(stranded_piles, consensus,
                                         mutant_consensus, mutant_samples)

  cat("Determining unique vs. shared variants \n")
  unique_sites <- determine_unique(mut_cov_counts, seq_err)

  cat("Gathering sample metadata \n")
  mutant_sample_ids <- samples[mutant_samples]
  mutation_class <- get_mutation_class(mutant_consensus)

  cat("Calling variants \n")
  p <- ulapply(mut_cov_matrices, call_variant)

  cat("Tabulating data \n")
  test_table <- cbind(as.data.frame(species, stringsAsFactors = FALSE),
                      genotype, isolate, mutant_sample_ids, 1:bp, mean_coverage,
                      consensus,
                      mutant_consensus, mutation_class, max_mut_afs, diff_afs,
                      strand_bias, unique_sites, coverage_proportion, p,
                      stringsAsFactors = FALSE)
  colnames(test_table) <- c("species", "genotype", "isolate", "sample", "pos",
                           "coverage", "ref", "alt", "class", "af", "af_diff",
                           "strand_bias", "unique", "coverage_proportion",
                           "p_value")

  cat("Merging multinucleotide indels \n")
  test_table <- merge_significant_indels(test_table)
  test_table <- rename_small_deletions(test_table)
  write_output(test_table, dir = file.path(outdir, TEST_TABLE_DIR),
               genotype = genotype, type = "csv")
}

if (!interactive()) {
  opts <- docopt::docopt(doc)
  main(opts["PILEUP"], opts["seqerr"], opts["output_dir"])
}
