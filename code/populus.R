# Genotype and time of day shape the Populus drought response
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15242


library(Biobase)
library(GEOquery)
library(limma)
library(stringr)

gset <- getGEO("GSE15242", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL4359", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Preprocessing ----------------------------------------------------------------

eset <- gset

dim(eset)
plotDensities(eset, legend = FALSE)

# The data is log2 GC-RMA signal
sum(rowMeans(exprs(eset)) > 6)
plotDensities(eset[rowMeans(exprs(eset)) > 6, ], legend = FALSE)
eset <- eset[rowMeans(exprs(eset)) > 6, ]

pData(eset) <- pData(eset)[, c("title", "clone:ch1", "time of day:ch1", "treatment:ch1" )]
colnames(pData(eset)) <- c("title", "type", "time", "water")

# Only keep one of the 4 timepoints
eset <- eset[, pData(eset)[, "time"] == "mid-day"]

# Clean up names
pData(eset)[, "type"] <- tolower(pData(eset)[, "type"])
pData(eset)[, "water"] <- ifelse(pData(eset)[, "water"] == "well-watered",
                                 "normal", "drought")
pData(eset)[, "rep"] <- sprintf("r%d",
                                as.integer(str_sub(pData(eset)[, "title"], -1, -1)))
pData(eset)[, c("title", "time")] <- NULL
colnames(eset) <- paste(pData(eset)[, "type"],
                        pData(eset)[, "water"],
                        pData(eset)[, "rep"], sep = "_")
head(pData(eset))

saveRDS(eset, file = "../data/populus-eset.rds")

# Analysis ---------------------------------------------------------------------

if (!exists("eset")) {
  eset <- readRDS("../data/populus-eset.rds")
}

dim(eset)
table(pData(eset)[, c("type", "water")])

# Create single variable
group <- with(pData(eset), paste(type, water, sep = "."))
group <- factor(group)

# Create design matrix with no intercept
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Count the number of samples modeled by each coefficient
colSums(design)

# Load package
library(limma)

# Create a contrasts matrix
cm <- makeContrasts(type_normal = nm6.normal - dn34.normal,
                    type_drought = nm6.drought - dn34.drought,
                    water_nm6 = nm6.drought - nm6.normal,
                    water_dn34 = dn34.drought - dn34.normal,
                    interaction = (nm6.drought - nm6.normal) - (dn34.drought - dn34.normal),
                    levels = design)

# View the contrasts matrix
cm

# Load package
library(limma)

# Fit the model
fit <- lmFit(eset, design)

# Fit the contrasts
fit2 <- contrasts.fit(fit, contrasts = cm)

# Calculate the t-statistics for the contrasts
fit2 <- eBayes(fit2)

# Summarize results
results <- decideTests(fit2)
summary(results)
