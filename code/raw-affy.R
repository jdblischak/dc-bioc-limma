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

# The Populus data provides raw CEL files generated from the platform:
#
# GPL198 	[ATH1-121501] Affymetrix Arabidopsis ATH1 Genome Array
#
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE53nnn/GSE53990/suppl/GSE53990_RAW.tar

getGEOSuppFiles(GEO = "GSE53990")
untar("GSE53990/GSE53990_RAW.tar", exdir = "GSE53990")

affybatch <- ReadAffy(filenames = Sys.glob("GSE53990/*CEL.gz"), compress = TRUE)

# Take the average of all the probes for a given gene. Do not do any other
# pre-processing.
eset <- expresso(affybatch,
                 bg.correct = FALSE,
                 normalize = FALSE,
                 pmcorrect.method = "pmonly",
                 summary.method = "avgdiff")

x <- exprs(eset)

dim(x)

plotDensities(x, legend = FALSE)

x <- log(x)

plotDensities(x, legend = FALSE)

x <- normalizeBetweenArrays(x)

plotDensities(x, legend = FALSE)

x <- x[rowMeans(x) > 6, ]

plotDensities(x, legend = FALSE)

dim(x)
