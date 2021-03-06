% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_main.R
\name{AlleleProfileR.qc.simple}
\alias{AlleleProfileR.qc.simple}
\title{Basic QC function for BAM data.frames}
\usage{
AlleleProfileR.qc.simple(bamtable)
}
\arguments{
\item{bamtable}{Data.frame containing the reads, as generated by AlleleProfileR.bam.}
}
\value{
bam data.frame
}
\description{
This function will discard low qual reads based on the flags. However, note that this function only employs really basic filtering. More advanced read filtering or preprocessing may be required and can be accomplished using AlleleProfileR.preprocess and external software.
}
\author{
Arne Bruyneel
}
