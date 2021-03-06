#!/usr/bin/env Rscript

# Download GSE40289
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40289

library("dplyr")
# Need to load Biobase after dplyr b/c it exports rlang::exprs, which masks
# Biobase::exprs
library("Biobase")
library("GEOquery")
library("stringr")

# Download data ----------------------------------------------------------------

if (.Platform$OS.type != "windows") {
  options('download.file.method.GEOquery' = 'wget')
}

geo <- getGEO("GSE40289", GSEMatrix = TRUE, getGPL = TRUE)[[1]]

# Some of the probes are duplicated! I'll need to fix this to include more
# interpretable rownames
nrow(fData(geo)) == length(unique(fData(geo)$NAME))

# Create nice phenotype data columns -------------------------------------------

pheno <- pData(geo) %>%
  select(title) %>%
  mutate(genotype = ifelse(grepl("WT", title), "wt", "top2b"),
         treatment = ifelse(grepl("PBS", title), "pbs", "dox"),
         timepoint = ifelse(grepl("16hr", title), "16hr", "72hr"),
         rep = str_extract(title, "\\d$"),
         rep = paste0("r", rep))

samples <- paste(pheno$genotype, pheno$treatment, pheno$timepoint,
                 pheno$rep, sep = ".")

rownames(pheno) <- samples

varmetadata <- data.frame(labelDescription = c(
  "Full description",
  "Genotype: WT = wild type, Tob2b = Tob2b null",
  "Treatment: PBS (control) vs Doxorubicin",
  "Timepoint: 16 or 72 hours post treatment",
  "Replicate: r1, r2, r3"
))

pdata <- new("AnnotatedDataFrame", data = pheno, varMetadata = varmetadata)

# Create nice feature data columns ---------------------------------------------

feature <- fData(geo) %>%
  mutate(chrom = str_split_fixed(CHROMOSOMAL_LOCATION, ":", 2)[, 1],
         pos = str_split_fixed(CHROMOSOMAL_LOCATION, ":", 2)[, 2],
         start = as.numeric(str_split_fixed(pos, "-", 2)[, 1]),
         end = as.numeric(str_split_fixed(pos, "-", 2)[, 2])) %>%
  select(probe = NAME,
         chrom,
         start,
         end,
         symbol = GENE_SYMBOL,
         ensembl = ENSEMBL_ID,
         refseq = GB_ACC,
         entrez = GENE_ID,
         name = GENE_NAME)

fvarmetadata <- data.frame(LabelDescription = c(
  "Microarray probe ID",
  "Chromosome",
  "Start position of probe",
  "End position of probe",
  "Gene symbol",
  "Ensembl ID",
  "RefSeq ID",
  "Entrez ID",
  "Full name"
))

fdata <- new("AnnotatedDataFrame", data = feature, varMetadata = fvarmetadata)

# Add column names to assay data -----------------------------------------------

assay <- exprs(geo)
rownames(assay) <- 1:nrow(assay)
colnames(assay) <- samples

# Create ExpressionSet ---------------------------------------------------------

eset <- ExpressionSet(assayData = assay, phenoData = pdata, featureData = fdata)
# Remove 72hr timepoint
eset <- eset[, eset$timepoint == "16hr"]
colnames(eset) <- str_replace(colnames(eset), "\\.16hr", "")
pData(eset)$timepoint <- NULL

# Filter ExpressionSet ---------------------------------------------------------

# The ExpressionSet RDS is over DataCamp's 5 MB limit

# Remove any probes without an entrez ID
eset <- eset[!is.na(fData(eset)[, "entrez"]), ]

# Remove duplicated probes. How is this even possible?!!
duplicate_probes <- duplicated(fData(eset)[, "probe"])
sum(duplicate_probes)
eset <- eset[!duplicate_probes, ]

# Use probes as rownames now that they are unique
rownames(eset) <- fData(eset)[, "probe"]

# Minimize feature data
fData(eset) <- fData(eset)[, c("symbol", "entrez", "chrom")]

# Minimize phenotype data
pData(eset)[, "title"] <- NULL

saveRDS(eset, file = "../data/dox.rds")
