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
\section{Accessing the Data}{
This resource can be accessed via ExperimentHub:

\if{html}{\out{<div class="sourceCode">}}\preformatted{
library(ExperimentHub)
eh <- ExperimentHub()

# Query and load
myfiles <- query(eh, "DMRsegaldata")
array_type <- myfiles[["EH####"]]  # Replace #### with actual EH ID

# Or use helper function
array_type <- ExperimentHub::loadResource(eh, "DMRsegaldata", "array_type.rda")
}\if{html}{\out{</div>}}

Alternatively, for local/development use:

\if{html}{\out{<div class="sourceCode">}}\preformatted{
data(array_type)
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
array_type <- ExperimentHub::loadResource(eh, "DMRsegaldata", "array_type.rda")
}

# Via LazyData
data(array_type)
array_type
}
\keyword{datasets}
