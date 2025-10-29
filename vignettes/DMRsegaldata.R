## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.width = 7,
    fig.height = 5
)

## ----eval=FALSE---------------------------------------------------------------
# # Install from Bioconductor (once submitted)
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#     install.packages("BiocManager")
# }
# BiocManager::install("DMRsegaldata")

## ----load-package-------------------------------------------------------------
library(DMRsegaldata)

## ----beta-data----------------------------------------------------------------
# Load the beta values matrix
data(beta)

# Examine the structure
dim(beta)
class(beta)

# Preview the first few rows and columns
beta[1:5, 1:5]

# Summary statistics
summary(beta[1:100, 1])

## ----pheno-data---------------------------------------------------------------
# Load phenotype data
data(pheno)

# View the structure
str(pheno)
head(pheno)

# Summary of sample groups
table(pheno$Sample_Group)

# Summary of gender distribution
table(pheno$Gender)
table(pheno$Sample_Group, pheno$Gender)

# Age distribution
summary(pheno$Age)

## ----dmps-data----------------------------------------------------------------
# Load DMP results
data(dmps)

# Examine the structure
head(dmps)
dim(dmps)

# Summary of significance levels
summary(dmps$pval)
summary(dmps$pval_adj)

# Top 10 most significant DMPs
head(dmps, 10)

# Count of significant DMPs at different thresholds
sum(dmps$pval_adj < 0.05)
sum(dmps$pval_adj < 0.01)
sum(dmps$pval_adj < 0.001)

## ----array-type---------------------------------------------------------------
# Load array type information
data(array_type)
array_type

## ----exploration--------------------------------------------------------------
# Check sample consistency between beta and pheno
all(colnames(beta) == rownames(pheno))

# Calculate mean methylation per sample
mean_methylation <- colMeans(beta, na.rm = TRUE)
names(mean_methylation) <- colnames(beta)

# Compare mean methylation between groups
cancer_samples <- rownames(pheno)[pheno$Sample_Group == "cancer"]
normal_samples <- rownames(pheno)[pheno$Sample_Group == "normal"]

mean_cancer <- mean(mean_methylation[cancer_samples])
mean_normal <- mean(mean_methylation[normal_samples])

cat("Mean methylation in cancer samples:", round(mean_cancer, 4), "\n")
cat("Mean methylation in normal samples:", round(mean_normal, 4), "\n")

## ----visualization, fig.cap="Distribution of mean methylation by sample group"----
# Create a data frame for plotting
plot_data <- data.frame(
    Sample = colnames(beta),
    MeanMethylation = mean_methylation,
    Group = pheno[colnames(beta), "Sample_Group"],
    Age = pheno[colnames(beta), "Age"],
    Gender = pheno[colnames(beta), "Gender"]
)

# Boxplot comparing groups
boxplot(MeanMethylation ~ Group,
    data = plot_data,
    main = "Mean Methylation by Sample Group",
    xlab = "Sample Group", ylab = "Mean Beta Value",
    col = c("cancer" = "lightcoral", "normal" = "lightblue")
)

# Add individual points
stripchart(MeanMethylation ~ Group,
    data = plot_data,
    vertical = TRUE, method = "jitter",
    add = TRUE, pch = 19, col = "black"
)

## ----top-dmps, fig.cap="Methylation levels at top 3 DMPs"---------------------
# Get top 3 DMPs
top_dmps <- head(rownames(dmps), 3)

# Set up plotting area
par(mfrow = c(1, 3))

# Plot each DMP
for (cpg in top_dmps) {
    if (cpg %in% rownames(beta)) {
        cpg_values <- beta[cpg, ]
        boxplot(cpg_values ~ pheno$Sample_Group,
            main = cpg,
            xlab = "Group", ylab = "Beta Value",
            col = c("lightcoral", "lightblue")
        )
        stripchart(cpg_values ~ pheno$Sample_Group,
            vertical = TRUE, method = "jitter",
            add = TRUE, pch = 19, col = "black"
        )
    }
}

par(mfrow = c(1, 1))

