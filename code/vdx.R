# Analysis of breast cancer VDX data for videos

# Prepare data -----------------------------------------------------------------

library(breastCancerVDX)
data("vdx")
class(vdx)
dim(vdx)
head(pData(vdx))
head(fData(vdx))
x <- exprs(vdx)
f <- fData(vdx)
p <- pData(vdx)

f <- f[, c("Gene.symbol", "EntrezGene.ID", "Chromosome.location")]
colnames(f) <- c("symbol", "entrez", "chrom")

# Recode er as 0 = negative and 1 = positive
p[, "er"] <- ifelse(p[, "er"] == 0, "negative", "positive")
p <- p[, c("id", "age", "er")]

# Explore data -----------------------------------------------------------------

png("../figure/ch01/boxplot.png",
    width = 480, height = 480 * 0.75, units = "px")
boxplot(x[1, ] ~ p[, "er"], main = f[1, "symbol"])
dev.off()
eset <- ExpressionSet(assayData = x,
                      phenoData = AnnotatedDataFrame(p),
                      featureData = AnnotatedDataFrame(f))
dim(eset)
boxplot(exprs(eset)[1, ] ~ pData(eset)[, "er"],
        main = fData(eset)[1, "symbol"])

# limma pipeline ---------------------------------------------------------------

design <- model.matrix(~er, data = pData(eset))
head(design, 2)
colSums(design)
table(pData(eset)[, "er"])

library(limma)

fit <- lmFit(eset, design)
head(fit$coefficients, 3)
fit <- eBayes(fit)
head(fit$t, 3)
results <- decideTests(fit[, "erpositive"])
summary(results)

# group-means ------------------------------------------------------------------

design <- model.matrix(~0 + er, data = pData(eset))
head(design)
colSums(design)

library(limma)
cm <- makeContrasts(status = erpositive - ernegative,
                    levels = design)
cm

fit <- lmFit(eset, design)
head(fit$coefficients)
fit2 <- contrasts.fit(fit, contrasts = cm)
head(fit2$coefficients)
fit2 <- eBayes(fit2)
results <- decideTests(fit2)
summary(results)

topTable(fit2)

# Visualize results  -----------------------------------------------------------

stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")
dim(stats)

png("../figure/ch03/vis-slides-%03d.png")

hist(runif(10000))

hist(stats[, "P.Value"])

volcanoplot(fit2, highlight = 5, names = fit2$genes[, "symbol"])

dev.off()

# Enrichment -------------------------------------------------------------------

topTable(fit2, number = 3)

# 1000 genes (10% in gene set), 100 are DE (10% in gene set)
fisher.test(matrix(c(10, 100, 90, 900), nrow = 2))

# 1000 genes (10% in gene set), 100 are DE (30% in gene set)
fisher.test(matrix(c(30, 100, 70, 900), nrow = 2))

head(fit2$genes, 3)
entrez <- fit2$genes[, "entrez"]


enrich_kegg <- kegga(fit2, geneid = entrez, species = "Hs")
topKEGG(enrich_kegg, number = 4)

enrich_go <- goana(fit2, geneid = entrez, species = "Hs")
topGO(enrich_go, ontology = "BP", number = 3)
