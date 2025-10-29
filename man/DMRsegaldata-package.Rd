\name{DMRsegaldata-package}
\alias{DMRsegaldata-package}
\alias{DMRsegaldata}
\docType{package}
\title{
Example DNA Methylation Data for DMRsegal
}
\description{
The DMRsegaldata package provides example DNA methylation datasets for use with the DMRsegal package. It contains preprocessed Illumina HumanMethylation450K array data from 10 cancer and 10 normal/healthy samples, along with associated phenotype information, differentially methylated positions (DMPs), and array type metadata.

This data package is designed to support the DMRsegal vignette and examples, providing users with realistic methylation data to explore differential methylation analysis workflows.
}
\details{
\tabular{ll}{
Package: \tab DMRsegaldata\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
License: \tab GPL (>= 2)\cr
LazyData: \tab yes\cr
}

The package contains four main data objects:

\describe{
  \item{\code{beta}}{A matrix of DNA methylation beta values (proportion of methylation) for CpG sites across 20 samples. Rows represent CpG sites and columns represent samples (10 cancer: T1-T10, 10 normal: H1-H10).}
  
  \item{\code{pheno}}{A data frame containing phenotype information for the 20 samples, including sample group (cancer vs. normal), age, and gender.}
  
  \item{\code{dmps}}{A data frame of differentially methylated positions identified by limma analysis comparing cancer vs. normal samples, filtered at FDR < 0.05. Includes p-values, q-values, and adjusted p-values.}
  
  \item{\code{array_type}}{A character string indicating the Illumina array platform used ("450K").}
}

All datasets are provided as .rda files with xz compression and are loaded automatically via LazyData.
}
\section{Data Preparation}{
The methylation data has been preprocessed and quality controlled. The preprocessing included:
\itemize{
  \item Quality control of raw Illumina 450K array data
  \item Normalization of beta values
  \item Filtering of low-quality probes
  \item Calculation of differentially methylated positions using limma
}

The DMP analysis was performed using empirical Bayes moderated t-statistics on M-values (log2 ratio of methylated to unmethylated signal), with adjustment for age and gender covariates.
}
\section{Accessing Data}{
There are two methods to access data from this package:

\subsection{Method 1: Via ExperimentHub (Recommended)}{
Resources can be accessed through the standard ExperimentHub interface:

\if{html}{\out{<div class="sourceCode">}}\preformatted{
library(ExperimentHub)
eh <- ExperimentHub()

# Query for all DMRsegaldata resources
myfiles <- query(eh, "DMRsegaldata")
myfiles

# Load specific resources by name
beta <- myfiles[["EH####"]]      # Replace with actual EH ID
pheno <- myfiles[["EH####"]]     # Replace with actual EH ID

# Or use package helper functions
library(ExperimentHub)
eh <- ExperimentHub()

# List all resources in this package
ExperimentHub::listResources(eh, "DMRsegaldata")

# Load a specific resource
beta <- ExperimentHub::loadResource(eh, "DMRsegaldata", "beta.rda")
}\if{html}{\out{</div>}}
}

\subsection{Method 2: Via LazyData (Development)}{
During package development or when installed locally with data/ directory:

\if{html}{\out{<div class="sourceCode">}}\preformatted{
data(beta)
data(pheno)
data(dmps)
data(array_type)
}\if{html}{\out{</div>}}
}

For detailed examples and workflows, see the package vignette:

\code{browseVignettes("DMRsegaldata")}
}
\section{Integration with DMRsegal}{
This data package is designed to work seamlessly with the DMRsegal package for differential methylation region analysis. The data objects can be directly used as inputs to DMRsegal functions.

See the DMRsegal package documentation for analysis workflows using this data.
}
\author{
Vasileios Lemonidis

Maintainer: Vasileios Lemonidis <vasileios.lemonidis@uantwerpen.be>
}
\references{
\itemize{
  \item DMRsegal package: \url{https://github.com/CMG-UA/DMRsegal}
  \item Illumina HumanMethylation450K array documentation
  \item Bioconductor workflows for DNA methylation analysis
}
}
\seealso{
\code{\link{beta}}, \code{\link{pheno}}, \code{\link{dmps}}, \code{\link{array_type}}

Related packages:
\itemize{
  \item \code{DMRsegal} - Differential methylation region analysis
  \item \code{minfi} - Illumina methylation array analysis
  \item \code{limma} - Linear models for microarray and RNA-seq data
}
}
\examples{
\dontrun{
## Access via ExperimentHub (recommended for end users)
library(ExperimentHub)
eh <- ExperimentHub()

## Query for DMRsegaldata resources
myfiles <- query(eh, "DMRsegaldata")
myfiles

## List all resources in this package
ExperimentHub::listResources(eh, "DMRsegaldata")

## Load resources using helper functions
beta <- ExperimentHub::loadResource(eh, "DMRsegaldata", "beta.rda")
pheno <- ExperimentHub::loadResource(eh, "DMRsegaldata", "pheno.rda")
dmps <- ExperimentHub::loadResource(eh, "DMRsegaldata", "dmps.rda")
array_type <- ExperimentHub::loadResource(eh, "DMRsegaldata", "array_type.rda")
}

## Access via LazyData (for local/development use)
data(beta)
data(pheno)
data(dmps)
data(array_type)

## Explore the data
dim(beta)
head(beta[, 1:5])

## Examine phenotype data
head(pheno)
table(pheno$Sample_Group)

## View DMPs
head(dmps)
cat("Number of significant DMPs:", nrow(dmps), "\n")

## Check array type
array_type

## Calculate mean methylation per sample group
cancer_samples <- rownames(pheno)[pheno$Sample_Group == "cancer"]
normal_samples <- rownames(pheno)[pheno$Sample_Group == "normal"]

mean_cancer <- mean(colMeans(beta[, cancer_samples], na.rm = TRUE))
mean_normal <- mean(colMeans(beta[, normal_samples], na.rm = TRUE))

cat("Mean methylation (cancer):", round(mean_cancer, 4), "\n")
cat("Mean methylation (normal):", round(mean_normal, 4), "\n")

## Examine a top DMP
if (nrow(dmps) > 0) {
  top_cpg <- rownames(dmps)[1]
  if (top_cpg \%in\% rownames(beta)) {
    boxplot(beta[top_cpg, ] ~ pheno$Sample_Group,
            main = paste("Top DMP:", top_cpg),
            xlab = "Sample Group", ylab = "Beta Value",
            col = c("lightcoral", "lightblue"))
  }
}
}
\keyword{package}
\keyword{datasets}
