% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multinomial_variant_calling.R
\name{get_minor_allele_frequencies}
\alias{get_minor_allele_frequencies}
\title{Generate a table of minor allele frequencies for each genome position in a
sample}
\usage{
get_minor_allele_frequencies(pile, consensus)
}
\arguments{
\item{pile}{matrix; allele counts at each genome position}

\item{consensus}{The consensus sequence for the samples}
}
\value{
Numeric vector of allele frequencies, named by allele.
}
\description{
Generate a table of minor allele frequencies for each genome position in a
sample
}
