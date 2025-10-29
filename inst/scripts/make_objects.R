dir.create("../../data", showWarnings = FALSE)
beta <- as.matrix(read.csv("../extdata/beta.tsv.gz", row.names = 1, header = TRUE, sep = "\t"))
save(beta, file = "../../data/beta.rda", compress = "xz")
pheno <- read.csv("../extdata/pheno.tsv", row.names = 1, header = TRUE, sep = "\t")
save(pheno, file = "../../data/pheno.rda", compress = "xz")
if (!requireNamespace("limma", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("limma")
}
mvalues <- log2(beta / (1 - beta + 1e-6) + 1e-6)
covs <- pheno[, c("Age", "Gender")]
pheno <- pheno[, "Sample_Group"]
pheno_and_covs <- cbind(list(x = pheno), covs)
design <- model.matrix(~., data = pheno_and_covs)
suppressWarnings(suppressMessages({
    sink("/dev/null")
    not_collinear <- apply(!is.na(coef(limma::lmFit(mvalues, design))), 2, all)
    sink()
}))
design <- design[, not_collinear]

fit <- limma::lmFit(mvalues, design)
fit <- limma::eBayes(fit)

pcols <- paste("x", unique(pheno), sep = "")
pcols <- pcols[pcols %in% colnames(design)]
anova.res <- limma::topTable(fit, coef = pcols, number = Inf, sort.by = "none")
tab <- data.frame(
    intercept = anova.res$AveExpr,
    pval = anova.res$P.Value
)
tab$f <- anova.res$B
rownames(tab) <- rownames(anova.res)
p0 <- siggenes::pi0.est(tab$pval[!is.na(tab$pval)])$p0
tab$qval <- siggenes::qvalue.cal(tab$pval, p0)
o <- order(tab$pval)
dmps <- tab[o, , drop = FALSE]
dmps$pval_adj <- p.adjust(dmps$pval, method = "BH")
dmps <- dmps[dmps$pval_adj < 0.05, ]

save(dmps, file = "../../data/dmps.rda", compress = "xz")

array_type <- "450K"
save(array_type, file = "../../data/array_type.rda", compress = "xz")
