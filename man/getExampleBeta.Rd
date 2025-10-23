\name{getExampleBeta}
\alias{getExampleBeta}
\title{Path to example beta file}
\usage{getExampleBeta()}
\value{A character scalar with the full path to \code{beta.tsv.gz} in this package.}
\description{Return the path to the example beta matrix file bundled in \code{DMRsegaldata}.}
\examples{
path <- getExampleBeta()
file.exists(path)
}
