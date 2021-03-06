% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variant_filtering.R
\name{filter_variants}
\alias{filter_variants}
\title{Exclude nonsignificant, non-unique and high-strand-bias variants}
\usage{
filter_variants(fdr_table, max_strand_bias)
}
\arguments{
\item{fdr_table}{A test table with qvalues}

\item{max_strand_bias}{Numeric; Maximum strand bias allowed for included
variants}
}
\value{
A filtered variant table
}
\description{
Exclude nonsignificant, non-unique and high-strand-bias variants
}
