getExampleBeta <- function() {
  system.file("extdata", "beta.tsv.gz", package = "DMRSegaldata", mustWork = TRUE)
}

getExamplePheno <- function() {
  system.file("extdata", "pheno.tsv", package = "DMRSegaldata", mustWork = TRUE)
}

listDMRSegalDataFiles <- function() {
  dir(system.file("extdata", package = "DMRSegaldata", mustWork = TRUE), full.names = TRUE)
}
