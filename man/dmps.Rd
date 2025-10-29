\name{dmps}
\alias{dmps}
\docType{data}
\title{
Differentially Methylated Positions (DMPs)
}
\description{
A data frame containing differentially methylated positions (CpG sites) identified by comparing cancer and normal samples. The DMPs were identified using limma with adjusted p-value < 0.05.
}
\usage{data("dmps")}
\format{
  A data frame with differentially methylated CpG sites as rows and the following columns:
  \describe{
    \item{\code{intercept}}{average expression (M-value) across all samples}
    \item{\code{pval}}{raw p-value from the statistical test}
    \item{\code{f}}{F-statistic or B-statistic value}
    \item{\code{qval}}{q-value (FDR-adjusted p-value using Storey's method)}
    \item{\code{pval_adj}}{adjusted p-value using Benjamini-Hochberg method}
  }
}
\details{
The DMPs were identified using limma's empirical Bayes moderated t-statistics on M-values (log2 ratio of methylated to unmethylated signal). Only CpG sites with adjusted p-value (Benjamini-Hochberg) < 0.05 are included. The results are sorted by p-value.

\subsection{Analysis Details}{
The differential methylation analysis included:
\itemize{
  \item Conversion of beta values to M-values for statistical testing
  \item Linear modeling with adjustment for age and gender covariates
  \item Empirical Bayes moderation of t-statistics (limma)
  \item Multiple testing correction using Benjamini-Hochberg method
  \item Filtering at FDR < 0.05
}
}
}
\section{Accessing the Data}{
This resource can be accessed via ExperimentHub:

\if{html}{\out{<div class="sourceCode">}}\preformatted{
library(ExperimentHub)
eh <- ExperimentHub()

# Query and load
myfiles <- query(eh, "DMRsegaldata")
dmps <- myfiles[["EH####"]]  # Replace #### with actual EH ID

# Or use helper function
dmps <- ExperimentHub::loadResource(eh, "DMRsegaldata", "dmps.rda")
}\if{html}{\out{</div>}}

Alternatively, for local/development use:

\if{html}{\out{<div class="sourceCode">}}\preformatted{
data(dmps)
}\if{html}{\out{</div>}}
}
\source{
Derived from the beta matrix and pheno data using limma differential methylation analysis.
}
\examples{
\dontrun{
# Via ExperimentHub
library(ExperimentHub)
eh <- ExperimentHub()
dmps <- ExperimentHub::loadResource(eh, "DMRsegaldata", "dmps.rda")
}

# Via LazyData
data(dmps)
head(dmps)
nrow(dmps)

# Summary of significance
summary(dmps$pval_adj)
sum(dmps$pval_adj < 0.01)  # Very significant DMPs
sum(dmps$pval_adj < 0.001) # Highly significant DMPs
}
\keyword{datasets}
