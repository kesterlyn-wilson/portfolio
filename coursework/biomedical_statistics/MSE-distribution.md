---
title: "RBIF111_HW_Week7"
author: "Kesterlyn Wilson"
date: "`r Sys.Date()`"
output: html_document
---

```{r, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Question 1

Cross-validation tests how well a predictive model generalizes to new real data. It works by splitting a dataset into k roughly equal parts, or “folds.” The model trains on k – 1 folds and tests on the remaining fold, repeating the process until each fold has served once as the test set. The average performance across all folds gives an estimate of the model’s predictive accuracy.

For example, in 5-fold cross-validation, the data is divided into five subsets. The model trains on four subsets and tests on the fifth, cycling through all five combinations. This method uses all available data efficiently and reduces bias that can come from a single random train/test split (Brownlee, 2020).

To set it up in R, you can randomly assign each sample to a fold:

`{r}`

`set.seed(123)`

`k <- 5 fold_ids <- sample(rep(1:k, length.out = nrow(data)))`

`folds <- split(seq_len(nrow(data)), fold_ids)`

Then, each fold becomes a test set in turn while the remaining samples form the training set.

Potential issues include data leakage, where information from the test set accidentally influences the model during training, and class imbalance, where one group appears more often in certain folds. Another issue is overfitting, especially if the same data are used both to tune hyperparameters and to estimate final performance. Random variation between folds can also lead to instability if the dataset is small.

Despite these challenges, cross-validation remains one of the most reliable methods for estimating model performance because it balances bias and variance more effectively than a single hold-out test (Kuhn & Johnson, 2019).



## Question 2

```{r setup}
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
# I manually downloaded from GEO because it was too slow to download in R
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

***Adding this temporary log_count as a demo because log1p(counts) is way too large. This will only use the first 1000 genes rather that all 20,000.***

```{r}
# --- TEMPORARY: run smaller, faster prototype ---
log_counts <- log_counts[1:1000, ]
```

## Question 3

```{r cross-validation loop, echo=FALSE}

# Create 5 folds
set.seed(123)
k <- 5
n <- length(Y)
fold_ids <- sample(rep(1:k, length.out = n))
folds <- split(seq_len(n), fold_ids)


# Function to get best gene per fold
cv_best_gene <- function(log_counts, Y, folds) {
  all.best.xval.genes <- character(length(folds))
  mse_values <- numeric(length(folds))
  
  for (i in seq_along(folds)) {
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(Y), test_idx)
    
    train_x <- log_counts[, train_idx, drop = FALSE]
    train_y <- Y[train_idx]
    test_x  <- log_counts[, test_idx, drop = FALSE]
    test_y  <- Y[test_idx]
    
    # Function to get ANOVA p-value for each gene
    get_p <- function(x, y = train_y) {
      fit <- try(lm(y ~ x), silent = TRUE)
      if (inherits(fit, "try-error")) return(NA_real_)
      a <- anova(fit)
      return(a$`Pr(>F)`[1])
    }
    
    # Compute p-values for all genes
    pvals <- apply(train_x, 1, get_p)
    
    # Select best gene (lowest p-value)
    best_gene <- names(which.min(pvals))
    all.best.xval.genes[i] <- best_gene
    
    # Refit model on training data using the best gene
    x_train_best <- as.numeric(train_x[best_gene, ])
    x_test_best  <- as.numeric(test_x[best_gene, ])
    
    model <- lm(train_y ~ x_train_best)
    pred <- predict(model, newdata = data.frame(x_train_best = x_test_best))
    mse_values[i] <- mean((pred - test_y)^2)
  }
  
  return(list(best_genes = all.best.xval.genes, mse = mse_values))
}


```

***n_sim \<-- 250 is also too large and causes a fatal error on my PC. For visualization, I ran this with with n_sim set to 10 and 50.***

```{r bootstrap simulation, n=50}
set.seed(123)
n_sim <- 50
bootstrap_results <- vector("list", n_sim)

for (i in 1:n_sim) {
  cat("Running simulation", i, "\n")
  
  # Bootstrap resampling
  boot_idx <- sample(seq_along(Y), replace = TRUE)
  
  boot_counts <- log_counts[, boot_idx]
  boot_Y <- Y[boot_idx]
  
  # Create new 5 folds for this bootstrap
  k <- 5
  n <- length(boot_Y)
  fold_ids <- sample(rep(1:k, length.out = n))
  folds <- split(seq_len(n), fold_ids)
  
  # Run cross-validation
  bootstrap_results[[i]] <- cv_best_gene(boot_counts, boot_Y, folds)
}


```

```{r summary}
# Combine all best genes
all_genes <- unlist(lapply(bootstrap_results, function(x) x$best_genes))

# Frequency of each gene being selected
gene_freq <- sort(table(all_genes), decreasing = TRUE)

# Show top 10
head(gene_freq, 10)

```

```{r visualize MSE}
mse_values <- unlist(lapply(bootstrap_results, function(x) x$mse))

boxplot(mse_values,
        main = "MSE Distribution Across 50 Bootstrap Simulations",
        ylab = "Mean Squared Error")

```


