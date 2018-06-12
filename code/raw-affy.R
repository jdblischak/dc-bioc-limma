# Obtaining some raw data to practice pre-processing
# https://raw.githubusercontent.com/jdblischak/tb-suscept/32a1e676d90c15b2f56d650dbdc529f3257611ba/analysis/thuong2008.Rmd

# Steps:
#
# 1. log-transform
# 2. quantile normalize
# 3. filter

library(affy)
library(Biobase)
library(GEOquery)
library(limma)

# Arabidopsis ------------------------------------------------------------------

# The Arabidopsis data provides raw CEL files generated from the platform:
#
# GPL198 	[ATH1-121501] Affymetrix Arabidopsis ATH1 Genome Array
#
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53990/suppl/GSE53990_RAW.tar
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53990

rds <- "../data/arabidopsis-eset-raw.rds"
if (!file.exists(rds)) {
  download <- getGEOSuppFiles(GEO = "GSE53990", baseDir = tempdir())
  tarfile <- rownames(download)
  untar(tarfile, exdir = dirname(tarfile))

  affybatch <- ReadAffy(filenames = Sys.glob(file.path(dirname(tarfile), "*CEL.gz")),
                        compress = TRUE)

  # Take the average of all the probes for a given gene. Do not do any other
  # pre-processing.
  eset <- expresso(affybatch,
                   bg.correct = FALSE,
                   normalize = FALSE,
                   pmcorrect.method = "pmonly",
                   summary.method = "avgdiff")

  # Get the phenotype data from the processed data
  eset_processed <- readRDS("../data/arabidopsis-eset.rds")
  pData(eset) <- pData(eset_processed)

  saveRDS(eset, file = rds)
} else {
  eset <- readRDS(rds)
}

# Load package
library(limma)

# View the distribution of the raw data
plotDensities(eset, legend = FALSE)

# Log tranform
exprs(eset) <- log(exprs(eset))
plotDensities(eset, legend = FALSE)

# Quantile normalize
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
plotDensities(eset, legend = FALSE)

# Load package
library(limma)

# View the normalized gene expression levels
plotDensities(eset, legend = FALSE)
abline(v = 5)

# Determine the genes with mean expression level greater than 5
keep <- rowMeans(exprs(eset)) > 5
sum(keep)

# Filter the genes
eset <- eset[keep, ]
plotDensities(eset, legend = FALSE)

# Populus ----------------------------------------------------------------------

# The Populus data provides raw CEL files generated from the platform:
#
# GPL4359	[Poplar] Affymetrix Poplar Genome Array
#
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE15242&format=file
#
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15242

rds <- "../data/populus-eset-raw.rds"
if (!file.exists(rds)) {
  download <- getGEOSuppFiles(GEO = "GSE15242", baseDir = tempdir())
  tarfile <- rownames(download)
  untar(tarfile, exdir = dirname(tarfile))

  affybatch <- ReadAffy(filenames = Sys.glob(file.path(dirname(tarfile), "*CEL.gz")),
                        compress = TRUE)

  # Take the average of all the probes for a given gene. Do not do any other
  # pre-processing.
  eset <- expresso(affybatch,
                   bg.correct = FALSE,
                   normalize = FALSE,
                   pmcorrect.method = "pmonly",
                   summary.method = "avgdiff")

  # Get the phenotype data from the processed data
  eset_processed <- readRDS("../data/populus-eset.rds")
  pData(eset) <- pData(eset_processed)


  saveRDS(eset, file = rds)
} else {
  eset <- readRDS(rds)
}

# Load package
library(limma)

# View the distribution of the raw data
plotDensities(eset, legend = FALSE)

# Log tranform
exprs(eset) <- log(exprs(eset))
plotDensities(eset, legend = FALSE)

# Quantile normalize
exprs(eset) <- normalizeBetweenArrays(exprs(eset))
plotDensities(eset, legend = FALSE)

# Load package
library(limma)

# View the normalized gene expression levels
plotDensities(eset, legend = FALSE)
abline(v = 5)

# Determine the genes with mean expression level greater than 5
keep <- rowMeans(exprs(eset)) > 5
sum(keep)

# Filter the genes
eset <- eset[keep, ]
plotDensities(eset, legend = FALSE)

# Misc -------------------------------------------------------------------------

100 - 1
log(100) - log(1)
.1 - .001
log(.1) - log(.001)
