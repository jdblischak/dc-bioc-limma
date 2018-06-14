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

png("../figure/ch01/boxplot.png")
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
head(design, 3)
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

fisher.test(matrix(c(10, 10, 40, 40), nrow = 2))

fisher.test(matrix(c(30, 10, 20, 40), nrow = 2))

head(fit2$genes, 3)
entrez <- fit2$genes[, "entrez"]


enrich_kegg <- kegga(fit2, geneid = entrez, species = "Hs")
topKEGG(enrich_kegg, number = 4)

enrich_go <- goana(fit2, geneid = entrez, species = "Hs")
topGO(enrich_go, ontology = "BP", number = 3)
