\name{beta}
\alias{beta}
\docType{data}
\title{
DNA Methylation Beta Values Matrix
}
\description{
A matrix of DNA methylation beta values for CpG sites across 20 samples (10 healthy and 10 cancer samples). The matrix contains preprocessed beta values from Illumina 450K methylation array data.
}
\usage{data("beta")}
\format{
  A numeric matrix with CpG sites as rows and sample identifiers as columns. Beta values range from 0 to 1, representing the proportion of methylation at each CpG site.
}
\details{
The beta values have been preprocessed and quality controlled. The data includes measurements from 10 healthy samples (H1-H10) and 10 cancer samples (T1-T10).
}
\source{
Example data for the DMRsegal package vignette and examples.
}
\examples{
data(beta)
dim(beta)
head(beta[, 1:5])
}
\keyword{datasets}
