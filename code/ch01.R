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

# Improved simulation ----------------------------------------------------------

set.seed(12345)
create_exp_mat <- function(n1, n2, ng,
                           alpha_mean, beta_mean, epsilon_sd) {
  status <- c(rep(0, n1), rep(1, n2))
  ns <- length(status)
  status <- matrix(status, nrow = 1)

  alpha <- rnorm(ng, mean = alpha_mean, sd = 1)
  beta <- matrix(rnorm(ng, mean = beta_mean, sd = 1), ncol = 1)
  epsilon <- matrix(rnorm(ng * ns, mean = 0, sd = epsilon_sd),
                    nrow = ng, ncol = ns)
  Yg <- alpha + beta %*% status + epsilon
  return(Yg)
}

x <- rbind(
  # 30 non-DE genes with high variance
  create_exp_mat(n1 = 3, n2 = 3, ng = 30, alpha_mean = 10, beta_mean = -1:1, epsilon_sd = 3),
  # 30 non-DE genes with low variance
  create_exp_mat(n1 = 3, n2 = 3, ng = 30, alpha_mean = 10, beta_mean = -1:1, epsilon_sd = 1),
  # 10 upregulated DE genes with low variance
  create_exp_mat(n1 = 3, n2 = 3, ng = 10, alpha_mean = 10, beta_mean = 5, epsilon_sd = 1),
  # 10 upregulated DE genes with high variance
  create_exp_mat(n1 = 3, n2 = 3, ng = 10, alpha_mean = 10, beta_mean = 5, epsilon_sd = 3),
  # 10 downregulated DE genes with low variance
  create_exp_mat(n1 = 3, n2 = 3, ng = 10, alpha_mean = 10, beta_mean = -5, epsilon_sd = 1),
  # 10 downregulated DE genes with high variance
  create_exp_mat(n1 = 3, n2 = 3, ng = 10, alpha_mean = 10, beta_mean = -5, epsilon_sd = 3)
)
image(t(x))

n1 <- 3
n2 <- 3
status <- c(rep(0, n1), rep(1, n2))
# limma
library("limma")
design <- model.matrix(~status)
fit <- lmFit(x, design)
head(fit$coefficients)
fit <- eBayes(fit)
results <- decideTests(fit[, 2])
summary(results)
stats <- topTable(fit, coef = "status", number = nrow(fit), sort.by = "none")

# lm
p <- numeric(length = nrow(x))
for (i in 1:length(p)) {
  mod <- lm(x[i, ] ~ status)
  p[i] <- summary(mod)$coefficients[2, 4]
}

stats <- cbind(stats,
               sd = apply(x, 1, sd),
               lm = p.adjust(p, method = "BH"))
library("ggplot2")
ggplot(stats, aes(x = sd, y = logFC, color = adj.P.Val < 0.05)) +
  geom_point()
ggplot(stats, aes(x = sd, y = logFC, color = lm < 0.05)) +
  geom_point()
plot(stats$adj.P.Val, stats$lm)
table(limma = stats$adj.P.Val < 0.05, lm = stats$lm < 0.05)
which(stats$adj.P.Val < 0.05)
which(stats$lm < 0.05)

ggplot(stats, aes(x = logFC, y = -log10(P.Value), color = adj.P.Val < 0.05)) +
  geom_point()
ggplot(stats, aes(x = logFC, y = -log10(P.Value), color = lm < 0.05)) +
  geom_point()
