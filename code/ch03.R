#!/usr/bin/env Rscript

# Download GSE75816
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75816

library("dplyr")
library("edgeR")
library("stringr")

dge <- readDGE(files = Sys.glob("../data/counts/*"))

dge$samples <- dge$samples %>%
  mutate(id = basename(colnames(dge)),
         group = str_sub(id, start = 1, end = -2)) %>%
  select(id, group, lib.size, norm.factors)

colnames(dge) <- dge$samples$id

saveRDS(dge, file = "../data/ch03.rds")
