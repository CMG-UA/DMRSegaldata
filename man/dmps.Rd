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
}
\source{
Derived from the beta matrix and pheno data using limma differential methylation analysis.
}
\examples{
data(dmps)
head(dmps)
nrow(dmps)
}
\keyword{datasets}
