# Differential expression analysis for doxorubicin study

# Analayze microarray data from Zhang et al., 2012, which measured gene
# expression in hearts from wild type and Top2b null mice treated with
# doxorubicin or a control.

# Process features -------------------------------------------------------------

# Pre
library(Biobase)
eset <- readRDS("../data/dox.rds")
library(limma)

# Video

png("../figure/ch04/dox-densities.png")

# View the distribution of the raw data
plotDensities(eset, group = pData(eset)[, "genotype"], legend = "topright")

dev.off()

# Exercise

# Log tranform
exprs(eset) <- log(exprs(eset))
plotDensities(eset,  group = pData(eset)[, "genotype"], legend = "topright")

# Quantile normalize
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
plotDensities(eset,  group = pData(eset)[, "genotype"], legend = "topright")

# Determine the genes with mean expression level greater than 0
keep <- rowMeans(exprs(eset)) > 0
sum(keep)

# Filter the genes
eset <- eset[keep, ]
plotDensities(eset, group = pData(eset)[, "genotype"], legend = "topright")

# Boxplot of Top2b -------------------------------------------------------------

# Create a boxplot of the expression level of Top2b to confirm the null mice
# have lower levels of Top2b expression.

# Pre
library(Biobase)
eset <- readRDS("../data/dox.rds")
library(limma)
exprs(eset) <- log(exprs(eset))
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
keep <- rowMeans(exprs(eset)) > 0
eset <- eset[keep, ]

# Exercise

# Find the row which contains Top2b expression data
top2b <- which(fData(eset)[, "symbol"] == "Top2b")

# Plot Top2b expression versus genotype
boxplot(exprs(eset)[top2b, ] ~ pData(eset)[, "genotype"],
        main = fData(eset)[top2b, ])

# Check sources of variation ---------------------------------------------------

# Pre
library(Biobase)
eset <- readRDS("../data/dox.rds")
library(limma)
exprs(eset) <- log(exprs(eset))
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
keep <- rowMeans(exprs(eset)) > 0
eset <- eset[keep, ]

# Exercise

# Plot principal components labeled by genotype
plotMDS(eset, labels = pData(eset)[, "genotype"], gene.selection = "common")

# Plot principal components labeled by treatment
plotMDS(eset, labels = pData(eset)[, "treatment"], gene.selection = "common")

# Design matrix ----------------------------------------------------------------

# Pre
library(Biobase)
eset <- readRDS("../data/dox.rds")
library(limma)
exprs(eset) <- log(exprs(eset))
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
keep <- rowMeans(exprs(eset)) > 0
eset <- eset[keep, ]

# Exercise

# Create single variable
group <- with(pData(eset), paste(genotype, treatment, sep = "."))
group <- factor(group)

# Create design matrix with no intercept
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Count the number of samples modeled by each coefficient
colSums(design)

# Contrasts matrix -------------------------------------------------------------

# Pre
library(Biobase)
eset <- readRDS("../data/dox.rds")
library(limma)
exprs(eset) <- log(exprs(eset))
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
keep <- rowMeans(exprs(eset)) > 0
eset <- eset[keep, ]
group <- with(pData(eset), paste(genotype, treatment, sep = "."))
group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Exercise

# Create a contrasts matrix
cm <- makeContrasts(dox_wt = wt.dox - wt.pbs,
                    dox_top2b = top2b.dox - top2b.pbs,
                    interaction = (top2b.dox - top2b.pbs) - (wt.dox - wt.pbs),
                    levels = design)

# View the contrasts matrix
cm

# Test for differential expression ---------------------------------------------

# Pre
library(Biobase)
eset <- readRDS("../data/dox.rds")
library(limma)
exprs(eset) <- log(exprs(eset))
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
keep <- rowMeans(exprs(eset)) > 0
eset <- eset[keep, ]
group <- with(pData(eset), paste(genotype, treatment, sep = "."))
group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
cm <- makeContrasts(dox_wt = wt.dox - wt.pbs,
                    dox_top2b = top2b.dox - top2b.pbs,
                    interaction = (top2b.dox - top2b.pbs) - (wt.dox - wt.pbs),
                    levels = design)

# Exercise

# Fit the model
fit <- lmFit(eset, design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrasts = cm)

# Calculate the t-statistics for the contrasts
fit2 <- eBayes(fit2)

# Summarize results
results <- decideTests(fit2)
summary(results)

# Create a Venn diagram
vennDiagram(results)

# Video

png("../figure/ch04/dox-venn.png")
vennDiagram(results)
dev.off()

# Histogram of p-values --------------------------------------------------------

# Pre
library(Biobase)
eset <- readRDS("../data/dox.rds")
library(limma)
exprs(eset) <- log(exprs(eset))
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
keep <- rowMeans(exprs(eset)) > 0
eset <- eset[keep, ]
group <- with(pData(eset), paste(genotype, treatment, sep = "."))
group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
cm <- makeContrasts(dox_wt = wt.dox - wt.pbs,
                    dox_top2b = top2b.dox - top2b.pbs,
                    interaction = (top2b.dox - top2b.pbs) - (wt.dox - wt.pbs),
                    levels = design)
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(fit2)

# Exercise

# Obtain the summary statistics for the contrast dox_wt
stats_dox_wt <- topTable(fit2, coef = "dox_wt", number = nrow(fit2),
                         sort.by = "none")
# Obtain the summary statistics for the contrast dox_top2b
stats_dox_top2b <- topTable(fit2, coef = "dox_top2b", number = nrow(fit2),
                            sort.by = "none")
# Obtain the summary statistics for the contrast interaction
stats_interaction <- topTable(fit2, coef = "interaction", number = nrow(fit2),
                              sort.by = "none")

# Create histograms of the p-values for each contrast
hist(stats_dox_wt[, "P.Value"])
hist(stats_dox_top2b[, "P.Value"])
hist(stats_interaction[, "P.Value"])

# Volcano plot -----------------------------------------------------------------

# Pre
library(Biobase)
eset <- readRDS("../data/dox.rds")
library(limma)
exprs(eset) <- log(exprs(eset))
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
keep <- rowMeans(exprs(eset)) > 0
eset <- eset[keep, ]
group <- with(pData(eset), paste(genotype, treatment, sep = "."))
group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
cm <- makeContrasts(dox_wt = wt.dox - wt.pbs,
                    dox_top2b = top2b.dox - top2b.pbs,
                    interaction = (top2b.dox - top2b.pbs) - (wt.dox - wt.pbs),
                    levels = design)
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(fit2)

# Exercise

# Extract the gene symbols
gene_symbols <- fit2$genes[, "symbol"]

# Create a volcano plot for the contrast dox_wt
volcanoplot(fit2, coef = "dox_wt", highlight = 5, names = gene_symbols)

# Create a volcano plot for the contrast dox_top2b
volcanoplot(fit2, coef = "dox_top2b", highlight = 5, names = gene_symbols)

# Create a volcano plot for the contrast interaction
volcanoplot(fit2, coef = "interaction", highlight = 5, names = gene_symbols)

# Pathway enrichment -----------------------------------------------------------

# Pre
library(Biobase)
eset <- readRDS("../data/dox.rds")
library(limma)
exprs(eset) <- log(exprs(eset))
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
keep <- rowMeans(exprs(eset)) > 0
eset <- eset[keep, ]
group <- with(pData(eset), paste(genotype, treatment, sep = "."))
group <- factor(group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
cm <- makeContrasts(dox_wt = wt.dox - wt.pbs,
                    dox_top2b = top2b.dox - top2b.pbs,
                    interaction = (top2b.dox - top2b.pbs) - (wt.dox - wt.pbs),
                    levels = design)
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrasts = cm)
fit2 <- eBayes(fit2)

# Exercise

# Extract the entrez gene IDs
entrez <- fit2$genes[, "entrez"]

# Test for enriched KEGG Pathways for contrast dox_wt
enrich_dox_wt <- kegga(fit2, coef = "dox_wt", geneid = entrez, species = "Mm")

# View the top 5 enriched KEGG pathways
topKEGG(enrich_dox_wt, number = 5)

# Test for enriched KEGG Pathways for contrast interaction
enrich_interaction <- kegga(fit2, coef = "interaction", geneid = entrez, species = "Mm")

# View the top 5 enriched KEGG pathways
topKEGG(enrich_interaction, number = 5)

# Conclusion -------------------------------------------------------------------

library(Biobase)
eset <- readRDS("../data/dox.rds")
library(limma)

png("../figure/ch04/dox-densities-all.png", width = 480 * 4)
par(mfrow = c(1, 4))
plotDensities(eset, group = pData(eset)[, "genotype"], legend = "topright")
exprs(eset) <- log(exprs(eset))
plotDensities(eset,  group = pData(eset)[, "genotype"], legend = "topright")
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
plotDensities(eset,  group = pData(eset)[, "genotype"], legend = "topright")
keep <- rowMeans(exprs(eset)) > 0
eset <- eset[keep, ]
plotDensities(eset, group = pData(eset)[, "genotype"], legend = "topright")
dev.off()

png("../figure/ch04/boxplot-top2b.png")
top2b <- which(fData(eset)[, "symbol"] == "Top2b")
boxplot(exprs(eset)[top2b, ] ~ pData(eset)[, "genotype"],
        main = fData(eset)[top2b, ])
dev.off()

png("../figure/ch04/dox-mds.png", width = 480 * 2)
par(mfrow = c(1, 2))
# Plot principal components labeled by genotype
plotMDS(eset, labels = pData(eset)[, "genotype"], gene.selection = "common")

# Plot principal components labeled by treatment
plotMDS(eset, labels = pData(eset)[, "treatment"], gene.selection = "common")
dev.off()

png("../figure/ch04/dox-hist.png", width = 480 * 3)
par(mfrow = c(1, 3))
stats_dox_wt <- topTable(fit2, coef = "dox_wt", number = nrow(fit2),
                         sort.by = "none")
stats_dox_top2b <- topTable(fit2, coef = "dox_top2b", number = nrow(fit2),
                            sort.by = "none")
stats_interaction <- topTable(fit2, coef = "interaction", number = nrow(fit2),
                              sort.by = "none")
hist(stats_dox_wt[, "P.Value"], main = "dox_wt")
hist(stats_dox_top2b[, "P.Value"], main = "dox_top2b")
hist(stats_interaction[, "P.Value"], main = "interaction")
dev.off()

png("../figure/ch04/dox-volcano.png", width = 480 * 3)
par(mfrow = c(1, 3), cex = 1.25)
# Extract the gene symbols
gene_symbols <- fit2$genes[, "symbol"]

# Create a volcano plot for the contrast dox_wt
volcanoplot(fit2, coef = "dox_wt", highlight = 5, names = gene_symbols)

# Create a volcano plot for the contrast dox_top2b
volcanoplot(fit2, coef = "dox_top2b", highlight = 5, names = gene_symbols)

# Create a volcano plot for the contrast interaction
volcanoplot(fit2, coef = "interaction", highlight = 5, names = gene_symbols)
dev.off()
