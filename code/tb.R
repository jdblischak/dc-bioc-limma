library("dplyr")
library("edgeR")
library("limma")
library("Biobase")

# Post batch correction --------------------------------------------------------

x_raw_file <- "table-s1.txt"
if (!file.exists(x_raw_file)) {
  download.file(url = "https://bitbucket.org/jdblischak/tb-data/raw/bc0f64292eb8c42a372b3a2d50e3d871c70c202e/table-s1.txt",
                destfile = x_raw_file)
}
x_raw <- read.delim(x_raw_file, stringsAsFactors = FALSE)

dim(x_raw)
x <- x_raw %>% select(starts_with("M"))
dim(x)
rownames(x) <- x_raw$id

f <- x_raw[, "name", drop = FALSE]
rownames(f) <- x_raw$id

p <- read.delim("https://bitbucket.org/jdblischak/tb-data/raw/bc0f64292eb8c42a372b3a2d50e3d871c70c202e/annotation.txt",
                   stringsAsFactors = FALSE)
dim(p)
rownames(p) <- paste(p$ind, p$bact, p$time, sep = ".")
colnames(x) <- rownames(p)

eset <- ExpressionSet(assayData = as.matrix(x),
                      phenoData = AnnotatedDataFrame(p),
                      featureData = AnnotatedDataFrame(f))


eset <- eset[, pData(eset)[, "bact"] %in% c("none", "Rv", "Salm") &
               pData(eset)[, "time"] == 18]

plotDensities(eset, legend = FALSE)
plotMDS(eset, labels = pData(eset)[, "bact"], gene.selection = "common")
plotMDS(eset, labels = pData(eset)[, "extr"], gene.selection = "common")

# Pre batch correction ---------------------------------------------------------

# Code adapted from https://bitbucket.org/jdblischak/tb/src/da33186481188c173f4d96fe8f4286cc3bbc1b9a/batch-correction-comparison.Rmd?at=master&fileviewer=file-view-default#batch-correction-comparison.Rmd-51
full_file <- "counts_per_sample.txt"
if (!file.exists(full_file)) {
  download.file(url = "https://bitbucket.org/jdblischak/tb-data/raw/bc0f64292eb8c42a372b3a2d50e3d871c70c202e/counts_per_sample.txt",
                destfile = full_file)
}
full <- read.delim(full_file, stringsAsFactors = FALSE)
dim(full)
full <- full[order(full$dir), ]
rownames(full) <- paste(full$ind, full$bact, full$time, sep = ".")
counts <- t(full[, grep("ENSG", colnames(full))])
# Filter lowly expressed genes
counts <- counts[rowSums(cpm(counts) > 1) > 6, ]
dim(counts)
stopifnot(colnames(counts) == rownames(p))

eset <- ExpressionSet(assayData = counts,
                      phenoData = AnnotatedDataFrame(p),
                      featureData = AnnotatedDataFrame(f))

exprs(eset) <- cpm(eset, log = TRUE)

eset <- eset[, pData(eset)[, "bact"] %in% c("none", "Rv", "GC") &
               pData(eset)[, "time"] == 18]

plotDensities(eset, legend = FALSE)
plotMDS(eset, labels = pData(eset)[, "bact"], gene.selection = "common")
plotMDS(eset, labels = pData(eset)[, "extr"], gene.selection = "common")
plotMDS(eset, labels = pData(eset)[, "time"], gene.selection = "common")

pca <- prcomp(t(exprs(eset)), scale. = TRUE)
plot(pca$x[, 1], pca$x[, 2], col = as.factor(pData(eset)[, "bact"]))
plot(pca$x[, 1], pca$x[, 2], col = pData(eset)[, "time"])
plot(pca$x[, 1], pca$x[, 2], col = pData(eset)[, "extr"])

library("ggplot2")
theme_set(theme_classic())
d <- cbind(pca$x, pData(eset))
d$infection <- ifelse(d$bact == "none", "control", "infection")
d$infection <- factor(d$infection, levels = c("infection", "control"))
p <- ggplot(d, aes(x = PC1, y = PC2, color = as.factor(time), size = infection)) +
  geom_point()
p
p + facet_wrap(~time)
p %+% d[d$time == 18, ]
