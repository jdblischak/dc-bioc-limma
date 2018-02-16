#!/usr/bin/env Rscript

# Create some synthetic data for illustrating concepts in chapter 01.

set.seed(1)
n <- 6
g <- 50
group <- rep(c("con", "treat"), each = n / 2)
samples <- paste0(group, 1:3)
genes <- sprintf("gene%02d", 1:g)
gexp <- matrix(rnorm(n = n * g, mean = 0, sd = 3),
               nrow = g, ncol = n,
               dimnames = list(genes, samples))
# Add mean expression level
means <- runif(n = g, min = 14, max = 16)
gexp <- sweep(gexp, 1, means, FUN = "+")
# Make 10% of genes DE
genes_de <- sample(1:g, size = round(0.25 * g))
for (i in genes_de) {
  gexp[i, ] <- ifelse(group == "treat", gexp[i, ] * 2, gexp[i, ])
}
heatmap(gexp)

saveRDS(gexp, "../data/ch01.rds")
