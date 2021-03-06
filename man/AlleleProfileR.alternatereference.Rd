% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_main.R
\name{AlleleProfileR.alternatereference}
\alias{AlleleProfileR.alternatereference}
\title{Determine the allele name of HDR sequences}
\usage{
AlleleProfileR.alternatereference(config, alternate, overwrite = T)
}
\arguments{
\item{config}{Configuration list, as generated by AlleleProfileR.setup.}

\item{alternate}{Path to .fasta file with HDR sequences}

\item{overwrite}{overwrite existing alternate.bam file}
}
\value{
Data.frame with the HDR variant names for each gene in the gene table.
}
\description{
This function determines the name of the variant that would be introduced if homologous recombination with the template were to occur.
}
\author{
Arne Bruyneel
}
