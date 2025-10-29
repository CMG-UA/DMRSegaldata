## This script creates the metadata.csv file for ExperimentHub submission
## The metadata follows the format required by ExperimentHubData::makeExperimentHubMetadata

## Create metadata data frame
metadata <- data.frame(
  ## Title: Short descriptive title for each resource
  Title = c(
    "DNA Methylation Beta Values Matrix for DMRsegal",
    "Phenotype Data for DMRsegal Example",
    "Differentially Methylated Positions for DMRsegal",
    "Array Type for DMRsegal Example Data"
  ),
  
  ## Description: Longer description of the resource
  Description = c(
    "Matrix of DNA methylation beta values from Illumina 450K array for 10 cancer and 10 normal samples. Beta values range from 0 to 1 representing proportion of methylation at each CpG site.",
    "Data frame containing sample phenotype information including sample group (cancer vs normal), age, and gender for 20 samples.",
    "Data frame of differentially methylated CpG sites identified by limma analysis comparing cancer vs normal samples, filtered at FDR < 0.05.",
    "Character string indicating the Illumina array type (450K) used for the methylation data."
  ),
  
  ## BiocVersion: Bioconductor version when resource was added
  BiocVersion = c("3.20", "3.20", "3.20", "3.20"),
  
  ## Genome: Genome build
  Genome = c("hg19", "hg19", "hg19", "hg19"),
  
  ## SourceType: Format of original data (e.g., FASTA, BAM, BigWig, etc.)
  SourceType = c("RDA", "RDA", "RDA", "RDA"),

  SourceUrl = c(
    "https://zenodo.org/records/17475298/files/beta.rda",
    "https://zenodo.org/records/17475298/files/pheno.rda",
    "https://zenodo.org/records/17475298/files/dmps.rda",
    "https://zenodo.org/records/17475298/files/array_type.rda"
  ),

  RDataPath = c(
    "DMRsegaldata/beta.rda",
    "DMRsegaldata/pheno.rda",
    "DMRsegaldata/dmps.rda",
    "DMRsegaldata/array_type.rda"
  ),
  
  ## SourceVersion: Version or date of source data
  SourceVersion = c("Oct 29 2025", "Oct 29 2025", "Oct 29 2025", "Oct 29 2025"),
  
  ## Species: Species name
  Species = c("Homo sapiens", "Homo sapiens", "Homo sapiens", "Homo sapiens"),
  
  ## TaxonomyId: NCBI taxonomy ID
  TaxonomyId = c(9606, 9606, 9606, 9606),

  Coordinate_1_based = c("TRUE", "NA", "TRUE", "NA"),

  ## DataProvider: Name of entity providing the data
  DataProvider = c(
    "Center for Oncological Research, University of Antwerp",
    "Center for Oncological Research, University of Antwerp",
    "Center for Oncological Research, University of Antwerp",
    "Center for Oncological Research, University of Antwerp"
  ),
  
  ## Maintainer: Maintainer name and email
  Maintainer = c(
    "Vasileios Lemonidis <vasileios.lemonidis@uantwerpen.be>",
    "Vasileios Lemonidis <vasileios.lemonidis@uantwerpen.be>",
    "Vasileios Lemonidis <vasileios.lemonidis@uantwerpen.be>",
    "Vasileios Lemonidis <vasileios.lemonidis@uantwerpen.be>"
  ),
  
  ## RDataClass: Class of R object returned
  RDataClass = c("matrix", "data.frame", "data.frame", "character"),
  
  ## DispatchClass: How to load data (e.g., Rda, H5File, FilePath)
  DispatchClass = c("Rda", "Rda", "Rda", "Rda"),
  

  
  ## Tags: Searchable tags for the resource
  Tags = c(
    "DNAMethylation:Illumina450k:BetaValues:Cancer:ExperimentData",
    "DNAMethylation:Phenotype:Cancer:ExperimentData",
    "DNAMethylation:DifferentialMethylation:Cancer:ExperimentData",
    "DNAMethylation:Illumina450k:ArrayType:ExperimentData"
  ),
  
  stringsAsFactors = FALSE
)

## Write metadata to CSV file
write.csv(metadata, 
          file = "../../inst/extdata/metadata.csv",
          row.names = FALSE)

cat("Metadata file created successfully at inst/extdata/metadata.csv\n")
cat("Number of resources:", nrow(metadata), "\n")
print(metadata[, c("Title", "RDataClass", "DispatchClass")])
