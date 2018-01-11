#!/usr/bin/env Rscript

# Download GSE75816
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75816

library("edgeR")

dge <- readDGE(files = Sys.glob("../data/counts/*"))
samples <- basename(colnames(dge))
colnames(dge) <- samples
dge$samples$group <- samples
dge$samples$files <- NULL

saveRDS(dge, file = "../data/ch03.rds")
