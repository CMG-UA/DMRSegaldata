\name{getExamplePheno}
\alias{getExamplePheno}
\title{Path to example phenotype file}
\usage{getExamplePheno()}
\value{A character scalar with the full path to \code{pheno.tsv} in this package.}
\description{Return the path to the example phenotype file bundled in \code{DMRsegaldata}.}
\examples{
path <- getExamplePheno()
file.exists(path)
}
