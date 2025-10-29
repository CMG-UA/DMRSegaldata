\name{pheno}
\alias{pheno}
\docType{data}
\title{
Phenotype Data for DNA Methylation Samples
}
\description{
A data frame containing phenotype information for 20 samples, including sample group (cancer vs. normal), age, and gender.
}
\usage{data("pheno")}
\format{
  A data frame with 20 observations on the following 3 variables:
  \describe{
    \item{\code{Sample_Group}}{a character vector indicating sample group: "cancer" or "normal"}
    \item{\code{Age}}{a numeric vector containing the age of each individual}
    \item{\code{Gender}}{a character vector indicating gender: "m" (male) or "f" (female)}
  }
}
\details{
The phenotype data contains information for 10 cancer samples (T1-T10) and 10 normal/healthy samples (H1-H10). Row names correspond to sample identifiers.
}
\section{Accessing the Data}{
This resource can be accessed via ExperimentHub:

\if{html}{\out{<div class="sourceCode">}}\preformatted{
library(ExperimentHub)
eh <- ExperimentHub()

# Query and load
myfiles <- query(eh, "DMRsegaldata")
pheno <- myfiles[["EH####"]]  # Replace #### with actual EH ID

# Or use helper function
pheno <- ExperimentHub::loadResource(eh, "DMRsegaldata", "pheno.rda")
}\if{html}{\out{</div>}}

Alternatively, for local/development use:

\if{html}{\out{<div class="sourceCode">}}\preformatted{
data(pheno)
}\if{html}{\out{</div>}}
}
\source{
Example data for the DMRsegal package vignette and examples.
}
\examples{
\dontrun{
# Via ExperimentHub
library(ExperimentHub)
eh <- ExperimentHub()
pheno <- ExperimentHub::loadResource(eh, "DMRsegaldata", "pheno.rda")
}

# Via LazyData
data(pheno)
head(pheno)
table(pheno$Sample_Group)

# Check sample distribution
table(pheno$Sample_Group, pheno$Gender)
summary(pheno$Age)
}
\keyword{datasets}
