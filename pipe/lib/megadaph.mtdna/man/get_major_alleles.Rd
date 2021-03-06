% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multinomial_variant_calling.R
\name{get_major_alleles}
\alias{get_major_alleles}
\title{Get the major allele at each genome position for a single sample
WARNING: Returns the first value if tied}
\usage{
get_major_alleles(pile)
}
\arguments{
\item{pile}{matrix; allele counts at each genome position}
}
\value{
Character; Vector of major alleles, each allele one of
        {A, C, G, T, +, -}
}
\description{
Get the major allele at each genome position for a single sample
WARNING: Returns the first value if tied
}
