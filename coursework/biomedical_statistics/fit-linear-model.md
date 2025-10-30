---
title: "RBIF111 - Week 6"
author: "Kesterlyn Wilson"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Question 1

I fit linear models of age \~ gene expression for all genes in the GSE60424 dataset. The histogram of ANOVA p-values shows a left-skewed distribution, suggesting that some genes’ expression levels change with age. After Bonferroni correction, only a small subset remain significant, reflecting its conservative control of the family-wise error rate. The Benjamini–Hochberg (FDR) adjustment preserves more discoveries because it limits the expected proportion of false positives. Gene rankings among the most significant genes remained consistent across methods, but less-significant genes changed order after correction.

I fit a linear model for each gene using *age* as the continuous dependent variable and gene expression as the predictor. The ANOVA p-value for the slope term represents how strongly each gene’s expression varies with age.

The histogram of raw p-values showed a strong left-skew, indicating that most genes had high p-values (no relationship), while a small subset exhibited very low p-values, suggesting potential age-associated expression patterns. The ranked −log₁₀(p) plot confirmed this: only a few genes stood out sharply before the curve flattened, which is typical for high-dimensional gene expression data where only a minority of genes correlate with a continuous trait.

After applying the Bonferroni correction, nearly all points collapsed toward zero, leaving only a handful of genes with significant adjusted p-values. This is expected because Bonferroni multiplies each p-value by the number of tests (\~20,000), providing strict family-wise error control but often being overly conservative.

When I applied the Benjamini–Hochberg (FDR) correction, the overall shape remained similar, but a slightly larger group of genes retained significance. FDR correction is less stringent and instead controls the expected proportion of false positives among significant findings. The top-ranking genes stayed consistent across all three methods, while mid-significance genes shifted in order after correction.

Overall, these results suggest that a small number of genes show strong evidence of age-related expression changes, while the majority remain unaffected. The different correction methods highlight the trade-off between strict error control (Bonferroni) and greater sensitivity to true biological signals (FDR).

```{r}
# Setup
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)

# Download metadata (pheno)
gse <- getGEO("GSE60424", GSEMatrix = TRUE)
pheno <- pData(gse[[1]])
colnames(pheno)

# Download counts
counts <- read.delim("GSE60424_GEOSubmit_FC1to11_normalized_counts.txt.gz", row.names = 1, check.names = FALSE)


# In order to align metadata with count matrix, reset rownames of pheno to use the library IDs
rownames(pheno) <- pheno$title

# Now align
common <- intersect(colnames(counts), rownames(pheno))
counts <- counts[, common]
pheno  <- pheno[common, ]

# Check alignment
all(colnames(counts) == rownames(pheno))  # should return TRUE

# Extract age (convert to numeric safely)
pheno$age <- as.numeric(pheno$`age:ch1`)

# remove NAs from data
keep <- !is.na(pheno$age)
counts <- counts[, keep]
pheno  <- pheno[keep, ]
Y <- pheno$age

# Log-transform counts (to stabilize variance)
log_counts <- log1p(counts)
```

```{r}
# ------------------------------------------------------------
#  Fit a linear model (age ~ gene expression) per gene
# ------------------------------------------------------------

# build a function that fits a linear model for each gene
get_p <- function(x, y = Y) {
  # try to fit a linear model of age (y) vs. gene expression (x)
  fit <- try(lm(y ~ x), silent = TRUE)
  # If the model fails (e.g., due to missing data or constant values), return NA
  if (inherits(fit, "try-error")) return(NA_real_)
  # Run ANOVA on the fitted model to get the p-value for the slope term
  a <- anova(fit)
  # Extract the p-value for the "x" variable (1st row of ANOVA table)
  return(a$`Pr(>F)`[1])
}

# Apply the function to each gene (row of the log-transformed count matrix)
pvals <- apply(log_counts, 1, get_p)

# Remove any NA values that may have resulted from failed model fits
pvals <- pvals[!is.na(pvals)]

# Plot a histogram showing the overall distribution of p-values
hist(pvals, breaks = 50,
     main = "Distribution of gene-wise ANOVA p-values",  # plot title
     xlab = "p-value",                                   # x-axis label
     col = "lightblue", border = "white")                # colors for the bars

# Identify the most and least significant genes
best_gene  <- names(which.min(pvals))   # gene with smallest p-value
worst_gene <- names(which.max(pvals))   # gene with largest p-value

# Print the names and p-values of those genes
cat("Most significant gene:", best_gene,  "p =", min(pvals), "\n")
cat("Least significant gene:", worst_gene, "p =", max(pvals), "\n")


```

```{r}
# ------------------------------------------------------------
# Multiple testing corrections
# ------------------------------------------------------------

# Apply Bonferroni correction (controls family-wise error rate)
p_bonf <- p.adjust(pvals, method = "bonferroni")

# Apply Benjamini–Hochberg (FDR) correction (controls expected false discovery rate)
p_bh   <- p.adjust(pvals, method = "BH")

# Combine all results into one data frame
results_part1 <- data.frame(
  GeneID = names(pvals),  # gene names from the rownames of log_counts
  p_raw  = pvals,         # raw p-values
  p_bonf = p_bonf,        # Bonferroni-corrected p-values
  p_bh   = p_bh           # FDR-corrected p-values
)

# (Optional) Save the results table to a CSV file for later use
write.csv(results_part1, "Part1_linear_model_results.csv", row.names = FALSE)
```

```{r}

# ------------------------------------------------------------
# Scatter plots of −log10(p) vs gene rank
# ------------------------------------------------------------

# Load ggplot2 for visualization
library(ggplot2)

# Define a helper function to make scatter plots of -log10(p) by rank
plot_ranked <- function(df, pcol, title) {
  # Order genes by their p-values
  df <- df[order(df[[pcol]]), ]
  # Assign rank based on sorted order
  df$rank <- seq_len(nrow(df))
  # Compute -log10(p) for y-axis; avoid log(0) using .Machine$double.xmin
  df$mlog10p <- -log10(pmax(df[[pcol]], .Machine$double.xmin))
  # Generate the scatter plot
  ggplot(df, aes(rank, mlog10p)) +
    geom_point(size = 0.6, alpha = 0.7, color = "steelblue") +   # dots
    labs(title = title,                                          # plot title
         x = "Genes ranked by significance",                     # x-axis label
         y = expression(-log[10](p))) +                          # y-axis label
    theme_bw()                                                   # clean background
}

# Create three plots for raw, Bonferroni, and FDR-adjusted p-values
p1_raw  <- plot_ranked(results_part1, "p_raw",  "Raw p-values")
p1_bonf <- plot_ranked(results_part1, "p_bonf", "Bonferroni correction")
p1_bh   <- plot_ranked(results_part1, "p_bh",   "FDR (BH) correction")

# Display each plot in sequence
print(p1_raw)
print(p1_bonf)
print(p1_bh)

```

## Question 2

I compared gene expression between the *MS pretreatment* and *Healthy Control* groups using a Mann–Whitney U test for each gene. This non-parametric test was chosen because it does not assume normality, which makes it suitable for RNA-seq data.

I also calculated an ad-hoc fold change as the ratio of average expression in the MS pretreatment samples to that in healthy controls.

```{r}
# ------------------------------------------------------------
# Part 2: Mann-Whitney U test — MS pretreatment vs Healthy Control
# ------------------------------------------------------------

# 1. Define full group factor
group_full <- as.factor(pheno$`diseasestatus:ch1`)
table(group_full)

# 2. Subset to the two groups of interest
keep <- group_full %in% c("Healthy Control", "MS pretreatment")
group <- droplevels(group_full[keep])
table(group)

# 3. Subset expression matrix and phenotype accordingly
X <- log_counts[, keep]
X <- as.matrix(X)

# 4. Sanitize group level names for safety in code
levels(group) <- make.names(levels(group))
levels(group)

# 5. Identify which is control and treatment
control_level   <- "Healthy.Control"
treatment_level <- "MS.pretreatment"
cat("Control:", control_level, " Treatment:", treatment_level, "\n")

# ------------------------------------------------------------
# Compute ad-hoc fold change and Mann-Whitney U test per gene
# ------------------------------------------------------------
mwu_results <- apply(X, 1, function(g) {
  g_ctrl <- g[group == control_level]
  g_trt  <- g[group == treatment_level]

  mean_ctrl <- mean(g_ctrl, na.rm = TRUE)
  mean_trt  <- mean(g_trt,  na.rm = TRUE)

  # ad-hoc fold change (treatment / control)
  fc <- mean_trt / mean_ctrl

  # Mann-Whitney U test (non-parametric)
  p <- tryCatch(
    wilcox.test(g_trt, g_ctrl, alternative = "two.sided")$p.value,
    error = function(e) NA_real_
  )

  return(c(fc = fc, p = p))
})

# 6. Build results data frame
mwu_df <- as.data.frame(t(mwu_results))
mwu_df$GeneID <- rownames(mwu_df)
mwu_df <- na.omit(mwu_df)

# 7. Multiple-testing corrections
mwu_df$p_bonf <- p.adjust(mwu_df$p, method = "bonferroni")
mwu_df$p_bh   <- p.adjust(mwu_df$p, method = "BH")

# 8. Sort and save
mwu_df <- mwu_df[order(mwu_df$p), ]
write.csv(mwu_df, "Part2_MWU_MSvsControl_results.csv", row.names = FALSE)



```

### Question 3

In this part of the analysis, I used permutation testing to validate the linear model results from Part 1 by comparing observed and randomized associations between gene expression and age. The most and least significant genes from the earlier analysis were selected for this step: **ENSG00000256618** (most significant) and **ENSG00000174437** (least significant). For each gene, I refit the linear model 2,000 times after randomly shuffling the age values to simulate a null hypothesis where there is no true association between expression and age. This process generated null distributions of both slope coefficients and ANOVA sum of squares values, which I compared to the observed statistics from the original data.

For *ENSG00000174437*, the null distributions of slopes were centered near zero, and the observed slope (indicated by the red dashed line) fell well within this range. Its observed variance explained was also similar to that expected by chance. The empirical permutation p-values for both slope and sum of squares were 1, confirming that this gene has no meaningful association with age. In contrast, *ENSG00000256618* showed a strong deviation from the null distributions. The observed slope and variance explained both fell far into the extreme right tail of their respective histograms, with empirical p-values of approximately 0. This indicates that the relationship between this gene’s expression and age is highly unlikely to occur by chance.

Overall, the permutation test confirmed the validity of the linear model approach. The most significant gene maintained a strong signal even under randomization, while the least significant gene behaved as expected under the null hypothesis. This demonstrates that permutation testing provides a robust, assumption-free method for assessing significance and supports the conclusion that a small subset of genes in this dataset genuinely vary with age.

```{r}
# ------------------------------------------------------------
# Permutation-based significance test
# ------------------------------------------------------------

# 1. Identify the most and least significant genes from Part 1
gene_top <- results_part1$GeneID[which.min(results_part1$p_bh)]   # smallest FDR p
gene_bot <- results_part1$GeneID[which.max(results_part1$p_bh)]   # largest FDR p

cat("Most significant gene:", gene_top, "\n")
cat("Least significant gene:", gene_bot, "\n")

# 2. Define a permutation testing function
perm_test_gene <- function(gene_id, Y, X, nperm = 2000) {
  # Extract expression for selected gene
  x <- as.numeric(X[gene_id, ])

  # Fit actual model (age ~ expression)
  fit <- lm(Y ~ x)
  slope_obs <- coef(fit)[["x"]]
  a <- anova(fit)
  ss_obs <- a$`Sum Sq`[1]   # variance explained by gene

  # Initialize storage vectors
  slopes <- numeric(nperm)
  ss <- numeric(nperm)

  # 3. Permute the age labels to break any true relationship
  for (i in seq_len(nperm)) {
    Y_perm <- sample(Y)                # shuffle ages
    fit_perm <- lm(Y_perm ~ x)         # fit permuted model
    a_perm <- anova(fit_perm)
    slopes[i] <- coef(fit_perm)[["x"]] # slope under null
    ss[i]     <- a_perm$`Sum Sq`[1]    # Sum Sq under null
  }

  # 4. Empirical p-values (two-sided for slope)
  p_slope <- mean(abs(slopes) >= abs(slope_obs))
  p_ss    <- mean(ss >= ss_obs)

  # 5. Return everything
  list(
    slope_obs = slope_obs,
    ss_obs = ss_obs,
    slopes = slopes,
    ss = ss,
    p_slope = p_slope,
    p_ss = p_ss
  )
}

# ------------------------------------------------------------
# Run permutation tests for the top and bottom genes
# ------------------------------------------------------------
set.seed(123)
res_top <- perm_test_gene(gene_top, Y, log_counts, nperm = 2000)
res_bot <- perm_test_gene(gene_bot, Y, log_counts, nperm = 2000)

# Display empirical p-values
cat("Top gene empirical p (slope):", res_top$p_slope,
    " | ANOVA SumSq p:", res_top$p_ss, "\n")
cat("Bottom gene empirical p (slope):", res_bot$p_slope,
    " | ANOVA SumSq p:", res_bot$p_ss, "\n")

# ------------------------------------------------------------
# Plot permutation distributions
# ------------------------------------------------------------
library(ggplot2)

plot_perm_hist <- function(values, obs, title, xlab) {
  ggplot(data.frame(v = values), aes(v)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    geom_vline(xintercept = obs, color = "red", linetype = "dashed", size = 1) +
    labs(title = title, x = xlab, y = "Count") +
    theme_bw()
}

# Top gene plots
print(plot_perm_hist(res_top$slopes, res_top$slope_obs,
                     paste0("Permutation slopes: ", gene_top),
                     "Slope under null"))
print(plot_perm_hist(res_top$ss, res_top$ss_obs,
                     paste0("Permutation Sum Sq: ", gene_top),
                     "Sum of Squares under null"))

# Bottom gene plots
print(plot_perm_hist(res_bot$slopes, res_bot$slope_obs,
                     paste0("Permutation slopes: ", gene_bot),
                     "Slope under null"))
print(plot_perm_hist(res_bot$ss, res_bot$ss_obs,
                     paste0("Permutation Sum Sq: ", gene_bot),
                     "Sum of Squares under null"))

```

![![](images/clipboard-1499701174.png)](images/clipboard-1977775718.png)

![](images/clipboard-2556718815.png)

![](images/clipboard-2290978163.png)
