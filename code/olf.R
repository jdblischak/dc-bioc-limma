library(Biobase)
library(limma)
# BiocInstaller::biocLite("HarmanData")
library(HarmanData)

data(OLF)
olf.info
table(olf.info)

colnames(olf.info) <- tolower(colnames(olf.info))
olf.info$treatment <- paste0("t", olf.info$treatment)
olf.info$batch <- paste0("b", olf.info$batch)

eset <- ExpressionSet(assayData = as.matrix(olf.data),
                      phenoData = AnnotatedDataFrame(olf.info))

saveRDS(eset, "../data/olf.rds")

dim(eset)
table(pData(eset))

plotMDS(eset, labels = pData(eset)[, "treatment"], gene.selection = "common")
plotMDS(eset, labels = pData(eset)[, "batch"], gene.selection = "common")

exprs(eset) <- removeBatchEffect(eset, batch = pData(eset)[, "batch"])

plotMDS(eset, labels = pData(eset)[, "treatment"], gene.selection = "common")
plotMDS(eset, labels = pData(eset)[, "batch"], gene.selection = "common")
