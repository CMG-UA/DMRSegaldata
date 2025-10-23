\name{array_type}
\alias{array_type}
\docType{data}
\title{
Methylation Array Type
}
\description{
A character string indicating the type of Illumina methylation array used for the data in this package.
}
\usage{data("array_type")}
\format{
  A character string with value "450K", indicating the data was generated using the Illumina HumanMethylation450 BeadChip array.
}
\details{
This variable specifies that the methylation data in this package comes from the Illumina 450K array platform, which measures methylation at approximately 450,000 CpG sites across the human genome.
}
\source{
Example data for the DMRsegal package vignette and examples.
}
\examples{
data(array_type)
array_type
}
\keyword{datasets}
