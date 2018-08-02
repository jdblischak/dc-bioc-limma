# Visualize and remove batch effects from TB data
#
# https://bitbucket.org/jdblischak/tb/src

library(dplyr)
library(limma)
library(edgeR)
# Have to load Biobase after dplyr so that exprs function works
library(Biobase)

# Download data ----------------------------------------------------------------

full_file <- "../data/counts_per_sample.txt"
full_file_url <- "https://bitbucket.org/jdblischak/tb-data/raw/bc0f64292eb8c42a372b3a2d50e3d871c70c202e/counts_per_sample.txt"
if (!file.exists(full_file)) {
  download.file(url = full_file_url, destfile = full_file)
}

# Pre-process data -------------------------------------------------------------

# Code adapted from https://bitbucket.org/jdblischak/tb/src/da33186481188c173f4d96fe8f4286cc3bbc1b9a/batch-correction-comparison.Rmd?at=master&fileviewer=file-view-default#batch-correction-comparison.Rmd-51

full <- read.delim(full_file, stringsAsFactors = FALSE)
dim(full)
full <- full[order(full$dir), ]
rownames(full) <- paste(full$ind, full$bact, full$time, sep = ".")

x <- t(full[, grep("ENSG", colnames(full))])
p <- full %>% select(ind, bact, time, extr, rin)
stopifnot(colnames(x) == rownames(p))

eset <- ExpressionSet(assayData = x,
                      phenoData = AnnotatedDataFrame(p))

# Filter lowly expressed genes
keep <- rowSums(cpm(exprs(eset)) > 1) > 6
sum(keep)
eset <- eset[keep, ]
dim(eset)

norm_factors <- calcNormFactors(exprs(eset))
exprs(eset) <- cpm(exprs(eset), lib.size = colSums(exprs(eset)) * norm_factors,
                   log = TRUE)

pData(eset)[, "infection"] <- ifelse(pData(eset)[, "bact"] == "none",
                                     "con", "inf")

pData(eset)[, "time"] <- ifelse(pData(eset)[, "time"] == 4,
                                "early", "late")

pData(eset)[, "batch"] <- sprintf("b%02d", pData(eset)[, "extr"])

table(pData(eset)[, c("time", "batch")])

# Analyze data -----------------------------------------------------------------

plotDensities(eset, legend = FALSE)

# Ch3 L2 plotMDS/removeBatchEffect
png("../figure/ch03/tb-pca-%03d.png")

par(cex = 1.25, # make the text larger
    cex.lab = 1.2, # make the axis titles larger (1.2 * 1.25 = 1.5)
    mar = c(4, 4, 1, 1) + 0.1 # Decrease the bottom and top margins
)

plotMDS(eset, labels = pData(eset)[, "time"], gene.selection = "common")

exprs(eset) <- removeBatchEffect(eset, batch = pData(eset)[, "batch"],
                                 covariates = pData(eset)[, "rin"])

plotMDS(eset, labels = pData(eset)[, "time"], gene.selection = "common")

dev.off()
