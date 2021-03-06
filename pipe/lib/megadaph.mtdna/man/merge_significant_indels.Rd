% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multinomial_variant_calling.R
\name{merge_significant_indels}
\alias{merge_significant_indels}
\title{Merge multinucleotide indels (p < 0.05) into single indel events
Deletions are renamed to correct VCF style format, but insertions will need
to be renamed manually.}
\usage{
merge_significant_indels(test_table)
}
\arguments{
\item{test_table}{data.frame; A table of variant information, rows are genome
positions}
}
\value{
data.frame; test_table with indels merged
}
\description{
Merge multinucleotide indels (p < 0.05) into single indel events
Deletions are renamed to correct VCF style format, but insertions will need
to be renamed manually.
}
