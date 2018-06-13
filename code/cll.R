# Analysis of leukemia CLL data for exercises

# Prepare data -----------------------------------------------------------------

library(CLL)
data("sCLLex")
class(sCLLex)
dim(sCLLex)
head(pData(sCLLex))
head(fData(sCLLex))
boxplot(exprs(sCLLex)[1, ] ~ pData(sCLLex)$Disease)
x <- exprs(sCLLex)
f <- fData(sCLLex)
p <- pData(sCLLex)
library(hgu95av2.db)
columns(hgu95av2.db)
# select(x, keys, columns, keytype, ...) Couldn't get multiVals argument to
# work. Seems to be ignored. I tried "first", "CharacterList", and a custom
# function to select the last possible value. None made any difference.
hope <- select(hgu95av2.db, keys = rownames(x),
               columns = c("SYMBOL", "ENTREZID", "MAP"))
hope <- hope[!duplicated(hope$PROBEID), ]
stopifnot(hope$PROBEID == rownames(x))
colnames(hope) <- c("probe", "symbol", "entrez", "chrom")
rownames(hope) <- hope$probe
hope[["probe"]] <- NULL
f <- hope

# Explore data -----------------------------------------------------------------

boxplot(x[1, ] ~ p[, "Disease"], main = f[1, "symbol"])
library(Biobase)
eset <- ExpressionSet(assayData = x,
                      phenoData = AnnotatedDataFrame(p),
                      featureData = AnnotatedDataFrame(f))
save(list = c("x", "f", "p"), file = "../data/cll.RData")
saveRDS(eset, file = "../data/cll-eset.rds")
dim(eset)

# Subset to only include the 1000th gene (row) and the first 10 samples
eset_sub <- eset[1000, 1:10]

# Check the dimensions of the subset
dim(eset_sub)

# Create a boxplot of the first gene in eset_sub
boxplot(exprs(eset_sub)[1, ] ~ pData(eset_sub)[, "Disease"], main = fData(eset_sub)[1, "symbol"])

# limma pipeline ---------------------------------------------------------------

# Create design matrix for CLL study
design <- model.matrix(~Disease, data = pData(eset))

# Count the number of samples modeled by each coefficient
colSums(design)

# Load package
library(limma)

# Fit the model
fit <- lmFit(eset, design)

# Calculate the t-statistics
fit <- eBayes(fit)

# Summarize results
results <- decideTests(fit[, "Diseasestable"])
summary(results)

# group-means ------------------------------------------------------------------

# Create design matrix with no intercept
design <- model.matrix(~0 + Disease, data = pData(eset))

# Count the number of samples modeled by each coefficient
colSums(design)

# Load package
library(limma)

# Create a contrasts matrix
cm <- makeContrasts(status = Diseaseprogres. - Diseasestable,
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

# Visualize results ------------------------------------------------------------

stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")
hist(stats[, "P.Value"])

volcanoplot(fit2, highlight = 4, names = fit2$genes$symbol)

# Enrichment -------------------------------------------------------------------

entrez <- fit2$genes$entrez
enrich_kegg <- kegga(fit2, geneid = entrez, species = "Hs")
topKEGG(enrich_kegg)

enrich_go <- goana(fit2, geneid = entrez, species = "Hs")
topGO(enrich_go, ontology = "BP")
