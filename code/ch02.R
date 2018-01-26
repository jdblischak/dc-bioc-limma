#!/usr/bin/env Rscript

# Download DAPARdata
# https://www.bioconductor.org/packages/release/data/experiment/html/DAPARdata.html
# https://www.ncbi.nlm.nih.gov/pubmed/27605098

# source("http://www.bioconductor.org/biocLite.R")
# biocLite("DAPARdata")

library("Biobase")
library("DAPARdata")
library("dplyr")
library("stringr")

data("Exp1_R25_prot")

# Create nice phenotype data columns -------------------------------------------

pheno <- pData(Exp1_R25_prot) %>%
  mutate(spikein = Label,
         rep = paste0("r", Bio.Rep)) %>%
  select(spikein, rep)

samples <- paste("spikein", pheno$spikein, pheno$rep, sep = ".")

rownames(pheno) <- samples

varmetadata <- data.frame(labelDescription = c(
  "Concentration (fmol) of spike-in of 48 human proteins",
  "Technical replicate: r1, r2, r3"
))

pdata <- new("AnnotatedDataFrame", data = pheno, varMetadata = varmetadata)

# Create nice feature data columns ---------------------------------------------

feature <- fData(Exp1_R25_prot) %>%
  select(protein = Protein.IDs)

rownames(feature) <- 1:nrow(feature)

fvarmetadata <- data.frame(LabelDescription = c(
  "Identification of all the proteins that a peptide mapped to"
))

fdata <- new("AnnotatedDataFrame", data = feature, varMetadata = fvarmetadata)

# Add column names to assay data -----------------------------------------------

assay <- exprs(Exp1_R25_prot)
rownames(assay) <- 1:nrow(assay)
colnames(assay) <- samples

# Export ExpressionSet ---------------------------------------------------------

eset <- ExpressionSet(assayData = assay, phenoData = pdata, featureData = fdata)
# Remove features which have NAs
missing <- apply(exprs(eset), 1, function(x) any(is.na(x)))
eset <- eset[!missing, ]
saveRDS(eset, file = "../data/ch02.rds")
