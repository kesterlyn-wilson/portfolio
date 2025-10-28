###############################################################################
# RBIF108 - Homework 1
# Author: Kesterlyn Wilson
# Title: Differential Expression and Pathway Analysis of RNA-seq Datasets
# Datasets: GSE202220, GSE213323, GSE160611
# Tools: DESeq2, ComplexHeatmap, fgsea, clusterProfiler, ReactomePA
###############################################################################

# =============================================================================
# Setup and Package Installation
# =============================================================================

# Check if BiocManager is installed; if not, install it
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required Bioconductor packages (only if not already installed)
required_packages <- c("DESeq2", "ComplexHeatmap", "fgsea", "qusage",
                       "EnhancedVolcano", "clusterProfiler", "ReactomePA", 
                       "org.Hs.eg.db", "enrichplot", "DOSE")

# Check which packages are not installed
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# Install only the missing packages
if(length(new_packages) > 0) {
  BiocManager::install(new_packages)
}

# Load all required libraries with error handling
tryCatch({
  library(DESeq2)           # Differential expression analysis
  library(data.table)       # Fast data manipulation
  library(dplyr)            # Data wrangling
  library(ComplexHeatmap)   # Advanced heatmap visualization
  library(fgsea)            # Fast gene set enrichment analysis
  library(qusage)           # Gene set analysis
  library(EnhancedVolcano)  # Volcano plots
  library(clusterProfiler)  # Functional enrichment analysis
  library(ReactomePA)       # Reactome pathway analysis
  library(org.Hs.eg.db)     # Human gene annotation database
  library(AnnotationDbi)    # Annotation database interface
  library(ggplot2)          # Publication-quality graphics
  library(gridExtra)        # Arrange multiple plots
  library(stringr)          # String manipulation
  library(enrichplot)       # Visualization for enrichment results
  library(DOSE)             # Disease ontology
  library(circlize)         # Circular visualization (for heatmap colors)
}, error = function(e) {
  stop("Error loading libraries: ", e$message)
})

# Set working directory to home folder
homedir <- "~"

# Create output directories for organized results
output_dir <- file.path(homedir, "RNAseq_results")
plots_dir <- file.path(output_dir, "plots")
tables_dir <- file.path(output_dir, "tables")
qc_dir <- file.path(output_dir, "quality_control")

# Create all directories if they don't exist
for (dir in c(output_dir, plots_dir, tables_dir, qc_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

cat("=============================================================================\n")
cat("RNA-seq Differential Expression Analysis Pipeline\n")
cat("=============================================================================\n\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Function to generate QC plots for DESeq2 results
generate_qc_plots <- function(dds, dataset_name) {
  cat(paste("Generating QC plots for", dataset_name, "...\n"))
  
  # Open PDF device for saving plots
  pdf(file.path(qc_dir, paste0(dataset_name, "_QC_plots.pdf")), width = 12, height = 8)
  
  # 1. Dispersion plot - shows relationship between mean and variance
  plotDispEsts(dds, main = paste(dataset_name, "- Dispersion Estimates"))
  
  # 2. Library size barplot - check if samples have similar sequencing depth
  lib_sizes <- colSums(counts(dds))
  barplot(lib_sizes/1e6, names.arg = colnames(dds), las = 2, 
          main = paste(dataset_name, "- Library Sizes (millions of reads)"),
          ylab = "Millions of reads", col = rainbow(length(lib_sizes)))
  
  # 3. Gene detection rate - how many genes detected per sample
  gene_detection <- colSums(counts(dds) > 0)
  barplot(gene_detection, names.arg = colnames(dds), las = 2,
          main = paste(dataset_name, "- Number of Detected Genes per Sample"),
          ylab = "Number of genes", col = rainbow(length(gene_detection)))
  
  dev.off()
  cat(paste("QC plots saved to:", qc_dir, "\n"))
}

# Function to create PCA plot for sample clustering visualization
create_pca_plot <- function(dds, dataset_name) {
  cat(paste("Creating PCA plot for", dataset_name, "...\n"))
  
  tryCatch({
    # Apply variance stabilizing transformation for PCA
    vsd <- vst(dds, blind = FALSE)
    
    # Generate PCA plot
    p <- plotPCA(vsd, intgroup = "Group") +
      ggtitle(paste(dataset_name, "- PCA Plot")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    
    # Save plot
    ggsave(file.path(plots_dir, paste0(dataset_name, "_PCA_plot.pdf")), 
           plot = p, width = 8, height = 6)
    
    return(p)
  }, error = function(e) {
    cat(paste("Warning: Could not generate PCA plot:", e$message, "\n"))
    return(NULL)
  })
}

# Function to create MA plot (log fold change vs mean expression)
create_ma_plot <- function(res, dataset_name) {
  cat(paste("Creating MA plot for", dataset_name, "...\n"))
  
  pdf(file.path(plots_dir, paste0(dataset_name, "_MA_plot.pdf")), width = 8, height = 6)
  plotMA(res, main = paste(dataset_name, "- MA Plot"), ylim = c(-5, 5))
  dev.off()
}

# Function to create volcano plot
create_volcano_plot <- function(res, dataset_name, padj_cutoff = 0.05, fc_cutoff = 1) {
  cat(paste("Creating volcano plot for", dataset_name, "...\n"))
  
  tryCatch({
    # Create enhanced volcano plot
    p <- EnhancedVolcano(res,
                        lab = res$Symbol,
                        x = 'log2FoldChange',
                        y = 'padj',
                        title = paste(dataset_name, '- Volcano Plot'),
                        pCutoff = padj_cutoff,
                        FCcutoff = fc_cutoff,
                        pointSize = 2.0,
                        labSize = 3.0,
                        legendPosition = 'right',
                        legendLabSize = 10,
                        col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                        colAlpha = 0.7,
                        drawConnectors = TRUE,
                        widthConnectors = 0.5,
                        max.overlaps = 20)
    
    # Save plot
    ggsave(file.path(plots_dir, paste0(dataset_name, "_volcano_plot.pdf")), 
           plot = p, width = 10, height = 8)
    
    return(p)
  }, error = function(e) {
    cat(paste("Warning: Could not generate volcano plot:", e$message, "\n"))
    return(NULL)
  })
}

# Function to perform pathway enrichment analysis
perform_pathway_enrichment <- function(results_df, dataset_name, 
                                      padj_cutoff = 0.05, fc_cutoff = 1) {
  cat(paste("Performing pathway enrichment analysis for", dataset_name, "...\n"))
  
  # Filter for significant DEGs
  sig_genes <- results_df %>%
    filter(padj < padj_cutoff, abs(log2FoldChange) > fc_cutoff) %>%
    filter(!is.na(ENTREZID))
  
  if (nrow(sig_genes) == 0) {
    cat("Warning: No significant genes found for enrichment analysis\n")
    return(NULL)
  }
  
  # Get upregulated and downregulated genes separately
  up_genes <- sig_genes %>% filter(log2FoldChange > 0) %>% pull(ENTREZID)
  down_genes <- sig_genes %>% filter(log2FoldChange < 0) %>% pull(ENTREZID)
  all_sig_genes <- sig_genes$ENTREZID
  
  # Background: all genes tested
  background <- results_df %>% filter(!is.na(ENTREZID)) %>% pull(ENTREZID)
  
  enrichment_results <- list()
  
  # 1. GO Biological Process enrichment
  tryCatch({
    ego_bp <- enrichGO(gene = all_sig_genes,
                       universe = background,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)
    
    if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
      enrichment_results$GO_BP <- ego_bp
      
      # Save results table
      write.table(as.data.frame(ego_bp), 
                  file.path(tables_dir, paste0(dataset_name, "_GO_BP_enrichment.txt")),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create dotplot
      pdf(file.path(plots_dir, paste0(dataset_name, "_GO_BP_dotplot.pdf")), 
          width = 10, height = 8)
      print(dotplot(ego_bp, showCategory = 20, title = paste(dataset_name, "- GO BP Enrichment")))
      dev.off()
    }
  }, error = function(e) {
    cat(paste("GO BP enrichment failed:", e$message, "\n"))
  })
  
  # 2. KEGG pathway enrichment
  tryCatch({
    kegg <- enrichKEGG(gene = all_sig_genes,
                       organism = 'hsa',
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
    
    if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
      # Convert gene IDs to symbols for readability
      kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      enrichment_results$KEGG <- kegg
      
      # Save results
      write.table(as.data.frame(kegg),
                  file.path(tables_dir, paste0(dataset_name, "_KEGG_enrichment.txt")),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create dotplot
      pdf(file.path(plots_dir, paste0(dataset_name, "_KEGG_dotplot.pdf")),
          width = 10, height = 8)
      print(dotplot(kegg, showCategory = 20, title = paste(dataset_name, "- KEGG Enrichment")))
      dev.off()
    }
  }, error = function(e) {
    cat(paste("KEGG enrichment failed:", e$message, "\n"))
  })
  
  # 3. Reactome pathway enrichment
  tryCatch({
    reactome <- enrichPathway(gene = all_sig_genes,
                              organism = "human",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              readable = TRUE)
    
    if (!is.null(reactome) && nrow(as.data.frame(reactome)) > 0) {
      enrichment_results$Reactome <- reactome
      
      # Save results
      write.table(as.data.frame(reactome),
                  file.path(tables_dir, paste0(dataset_name, "_Reactome_enrichment.txt")),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create dotplot
      pdf(file.path(plots_dir, paste0(dataset_name, "_Reactome_dotplot.pdf")),
          width = 10, height = 8)
      print(dotplot(reactome, showCategory = 20, title = paste(dataset_name, "- Reactome Enrichment")))
      dev.off()
    }
  }, error = function(e) {
    cat(paste("Reactome enrichment failed:", e$message, "\n"))
  })
  
  return(enrichment_results)
}

# Function to perform GSEA (Gene Set Enrichment Analysis)
perform_gsea <- function(results_df, dataset_name) {
  cat(paste("Performing GSEA for", dataset_name, "...\n"))
  
  tryCatch({
    # Prepare ranked gene list (rank by log2FC * -log10(pvalue))
    gene_list <- results_df %>%
      filter(!is.na(ENTREZID), !is.na(log2FoldChange), !is.na(pvalue)) %>%
      mutate(rank = log2FoldChange * -log10(pvalue)) %>%
      arrange(desc(rank))
    
    # Create named vector for GSEA
    ranks <- gene_list$rank
    names(ranks) <- gene_list$ENTREZID
    
    # Run GSEA with KEGG pathways
    gsea_kegg <- gseKEGG(geneList = ranks,
                         organism = 'hsa',
                         pvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         verbose = FALSE)
    
    if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
      # Make readable
      gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      
      # Save results
      write.table(as.data.frame(gsea_kegg),
                  file.path(tables_dir, paste0(dataset_name, "_GSEA_KEGG.txt")),
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Create enrichment plot for top pathways
      pdf(file.path(plots_dir, paste0(dataset_name, "_GSEA_enrichment_plot.pdf")),
          width = 12, height = 8)
      print(dotplot(gsea_kegg, showCategory = 15, title = paste(dataset_name, "- GSEA KEGG")))
      dev.off()
    }
    
    return(gsea_kegg)
  }, error = function(e) {
    cat(paste("GSEA failed:", e$message, "\n"))
    return(NULL)
  })
}

# =============================================================================
# SECTION 1: Differential Expression (DESeq2) – GSE202220
# =============================================================================

cat("\n=============================================================================\n")
cat("ANALYSIS 1: GSE202220\n")
cat("=============================================================================\n\n")

# Construct URL for downloading raw count data from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"

# Build full path to count data file
path <- paste(urld, "acc=GSE202220", 
              "file=GSE202220_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&")

# Read count matrix and convert to matrix format with GeneID as rownames
tbl <- as.matrix(data.table::fread(path, header = TRUE, colClasses = "integer"), 
                 rownames = "GeneID")

cat(paste("Downloaded count matrix with", nrow(tbl), "genes and", ncol(tbl), "samples\n"))

# Download and read gene annotation file (contains gene symbols and descriptions)
apath <- paste(urld, "type=rnaseq_counts", 
               "file=Human.GRCh38.p13.annot.tsv.gz", sep="&")
annot <- fread(apath, header = TRUE, quote = "", 
               stringsAsFactors = FALSE, data.table = FALSE)
rownames(annot) <- annot$GeneID

# Define experimental design
# gsms string: "1" = Treatment, "0" = Control
# Pattern "10101010" means alternating Tx and Ctrl samples
gsms <- "10101010"

# Convert string to factor for group assignment
gs <- factor(strsplit(gsms, split = "")[[1]])

# Create readable group names (make.names ensures valid R variable names)
groups <- make.names(c("Tx", "Ctrl"))

# Assign group labels to factor levels
levels(gs) <- groups

# Create sample information dataframe (required for DESeq2)
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

cat("Sample groups:\n")
print(table(sample_info$Group))

# Pre-filter low-count genes to improve statistical power
# Keep genes with at least 10 reads in minimum group size samples
keep <- rowSums(tbl >= 10) >= min(table(gs))
tbl <- tbl[keep, ]
cat(paste("Retained", nrow(tbl), "genes after filtering\n"))

# Create DESeq2 dataset object
# countData: raw count matrix
# colData: sample metadata
# design: formula specifying the experimental design
ds <- DESeqDataSetFromMatrix(countData = tbl, 
                              colData = sample_info, 
                              design = ~Group)

# Run DESeq2 differential expression analysis
# test="Wald": standard pairwise comparison
# sfType="poscount": normalization method suitable for datasets with extreme outliers
ds <- DESeq(ds, test = "Wald", sfType = "poscount")

# Generate QC plots
generate_qc_plots(ds, "GSE202220")

# Create PCA plot
create_pca_plot(ds, "GSE202220")

# Extract results for Tx vs Ctrl comparison
# alpha: significance threshold for independent filtering
# pAdjustMethod: multiple testing correction method
tT <- results(ds, contrast = c("Group", groups[1], groups[2]), 
              alpha = 0.05, pAdjustMethod = "fdr")

# Create MA plot
create_ma_plot(tT, "GSE202220")

# Merge DESeq2 results with gene annotations
tT <- merge(as.data.frame(tT), annot, by = 0, sort = FALSE)

# Select and reorder columns for final output
tT <- subset(tT, select = c("GeneID", "padj", "pvalue", "lfcSE", "stat",
                            "log2FoldChange", "baseMean", "Symbol", "Description"))

# Write results to tab-delimited file
write.table(tT, file = file.path(tables_dir, "GSE202220_DEGs.txt"), 
            row.names = FALSE, sep = "\t", quote = FALSE)

# Create volcano plot
create_volcano_plot(tT, "GSE202220")

# Prepare data for enrichment analysis
GSE202220_enrichment <- tT %>%
  mutate(Subset_Comparison = "GSE202220") %>%
  rename(ENTREZID = GeneID, SYMBOL = Symbol)

# Perform pathway enrichment
enrichment_202220 <- perform_pathway_enrichment(GSE202220_enrichment, "GSE202220")

# Perform GSEA
gsea_202220 <- perform_gsea(GSE202220_enrichment, "GSE202220")

cat("\nGSE202220 analysis complete!\n")

# =============================================================================
# SECTION 2: Differential Expression (DESeq2) – GSE213323
# =============================================================================

cat("\n=============================================================================\n")
cat("ANALYSIS 2: GSE213323\n")
cat("=============================================================================\n\n")

# Download GSE213323 count data
path <- paste(urld, "acc=GSE213323", 
              "file=GSE213323_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&")
tbl <- as.matrix(fread(path, header = TRUE, colClasses = "integer"), 
                 rownames = "GeneID")

cat(paste("Downloaded count matrix with", nrow(tbl), "genes and", ncol(tbl), "samples\n"))

# Use same annotation file as before
annot <- fread(apath, header = TRUE, quote = "", 
               stringsAsFactors = FALSE, data.table = FALSE)
rownames(annot) <- annot$GeneID

# Define experimental design for GSE213323
# "X" marks samples to exclude from analysis
# First 3 samples = Treatment, next 3 = Control
gsms <- "111000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

# Split string into individual characters
sml <- strsplit(gsms, split = "")[[1]]

# Select only non-"X" samples
sel <- which(sml != "X")
sml <- sml[sel]

# Subset count matrix to selected samples only
tbl <- tbl[, sel]

cat(paste("Selected", ncol(tbl), "samples for analysis\n"))

# Create factor for group assignment
gs <- factor(sml)
groups <- make.names(c("Tx", "Ctrl"))
levels(gs) <- groups

# Create sample metadata
sample_info <- data.frame(Group = gs, row.names = colnames(tbl))

cat("Sample groups:\n")
print(table(sample_info$Group))

# Filter low-count genes
keep <- rowSums(tbl >= 10) >= min(table(gs))
tbl <- tbl[keep, ]
cat(paste("Retained", nrow(tbl), "genes after filtering\n"))

# Create DESeq2 object
ds <- DESeqDataSetFromMatrix(countData = tbl, 
                              colData = sample_info, 
                              design = ~Group)

# Run differential expression analysis
ds <- DESeq(ds, test = "Wald", sfType = "poscount")

# Generate QC plots
generate_qc_plots(ds, "GSE213323")

# Create PCA plot
create_pca_plot(ds, "GSE213323")

# Extract and annotate results
tT <- results(ds, contrast = c("Group", groups[1], groups[2]), 
              alpha = 0.05, pAdjustMethod = "fdr")

# Create MA plot
create_ma_plot(tT, "GSE213323")

tT <- merge(as.data.frame(tT), annot, by = 0, sort = FALSE)
tT <- subset(tT, select = c("GeneID", "padj", "pvalue", "lfcSE", "stat",
                            "log2FoldChange", "baseMean", "Symbol", "Description"))

# Save results
write.table(tT, file = file.path(tables_dir, "GSE213323_DEGs.txt"), 
            row.names = FALSE, sep = "\t", quote = FALSE)

# Create volcano plot
create_volcano_plot(tT, "GSE213323")

# Prepare data for enrichment analysis
GSE213323_enrichment <- tT %>%
  mutate(Subset_Comparison = "GSE213323") %>%
  rename(ENTREZID = GeneID, SYMBOL = Symbol)

# Perform pathway enrichment
enrichment_213323 <- perform_pathway_enrichment(GSE213323_enrichment, "GSE213323")

# Perform GSEA
gsea_213323 <- perform_gsea(GSE213323_enrichment, "GSE213323")

cat("\nGSE213323 analysis complete!\n")

# =============================================================================
# SECTION 3: Merge Results from Both Datasets
# =============================================================================

cat("\n=============================================================================\n")
cat("MERGING RESULTS\n")
cat("=============================================================================\n\n")

# Combine both datasets into single table
DEseqrbound <- rbind(GSE202220_enrichment, GSE213323_enrichment)

# Write merged results to Excel-compatible format
fwrite(DEseqrbound, 
       file.path(tables_dir, "GSE202220_GSE213323_DEGAnalysis.xls"), 
       row.names = FALSE, quote = FALSE, sep = "\t")

cat("Merged results saved!\n")

# =============================================================================
# SECTION 4: Generate Comprehensive Summary Statistics
# =============================================================================

cat("\n=============================================================================\n")
cat("SUMMARY STATISTICS\n")
cat("=============================================================================\n\n")

# Count significant DEGs in each dataset (padj < 0.05)
sig_threshold <- 0.05
fc_threshold <- 1

summary_stats <- DEseqrbound %>%
  group_by(Subset_Comparison) %>%
  summarise(
    Total_Genes = n(),
    Significant_DEGs = sum(padj < sig_threshold, na.rm = TRUE),
    Highly_Significant = sum(padj < 0.01, na.rm = TRUE),
    Upregulated = sum(padj < sig_threshold & log2FoldChange > 0, na.rm = TRUE),
    Downregulated = sum(padj < sig_threshold & log2FoldChange < 0, na.rm = TRUE),
    Highly_Upregulated = sum(padj < sig_threshold & log2FoldChange > fc_threshold, na.rm = TRUE),
    Highly_Downregulated = sum(padj < sig_threshold & log2FoldChange < -fc_threshold, na.rm = TRUE),
    .groups = 'drop'
  )

# Print summary
cat("\nDifferential Expression Summary:\n")
print(summary_stats)

# Save summary statistics
write.table(summary_stats, 
            file = file.path(tables_dir, "DEG_Summary_Statistics.txt"),
            row.names = FALSE, sep = "\t", quote = FALSE)

# Create summary visualization
pdf(file.path(plots_dir, "DEG_summary_barplot.pdf"), width = 10, height = 6)
summary_long <- summary_stats %>%
  select(Subset_Comparison, Upregulated, Downregulated) %>%
  tidyr::pivot_longer(cols = c(Upregulated, Downregulated), 
                      names_to = "Direction", values_to = "Count")

p <- ggplot(summary_long, aes(x = Subset_Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_bw() +
  labs(title = "Differentially Expressed Genes by Dataset",
       x = "Dataset", y = "Number of DEGs") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
print(p)
dev.off()

# =============================================================================
# SECTION 5: Heatmap Function for Selected Pathways
# =============================================================================

# Function to create heatmap of genes from selected pathways
# This allows visualization of expression patterns across datasets
make_selected_pathway_heatmap <- function(selectAnnotate, 
                                          DEseqrbound,
                                          gene_col = "SYMBOL",
                                          comparison_col = "Subset_Comparison",
                                          value_col = "log2FoldChange",
                                          title = "Selected Pathway") {
  
  cat(paste("Creating heatmap for:", title, "\n"))
  
  tryCatch({
    # Parse gene IDs from pathway annotation (geneID field contains "/" separated values)
    spl <- str_split(selectAnnotate$geneID, "/")
    
    # Extract unique gene symbols, removing any fold change annotations
    selected <- unique(unlist(lapply(spl, function(x) {
      gsub("\\(logFC.+[0-9]\\)", "", x)
    })))
    
    # Filter DESeq results to only include genes in the selected pathway
    selDT <- as.data.table(DEseqrbound)[get(gene_col) %in% selected]
    
    if (nrow(selDT) == 0) {
      cat("Warning: No matching genes found for this pathway\n")
      return(NULL)
    }
    
    # Convert from long to wide format for heatmap
    # Each row = gene, each column = dataset comparison
    wide <- dcast(selDT, 
                  formula = as.formula(paste(gene_col, "~", comparison_col)),
                  value.var = value_col,
                  fun.aggregate = mean)  # Average if duplicate genes exist
    
    # Extract gene names for row labels
    gene_names <- wide[[gene_col]]
    
    # Remove gene name column and convert to matrix for plotting
    mat <- as.matrix(wide[, -1, with = FALSE])
    rownames(mat) <- gene_names
    
    # Create heatmap with ComplexHeatmap
    ht <- Heatmap(
      mat,
      name = "log2FC",  # Legend title
      column_title = title,  # Main title
      row_names_gp = gpar(fontsize = 8),  # Row label font size
      column_names_gp = gpar(fontsize = 10),  # Column label font size
      cluster_rows = TRUE,  # Cluster genes by similarity
      cluster_columns = FALSE,  # Don't cluster datasets
      show_row_names = TRUE,
      show_column_names = TRUE,
      col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),  # Color scale
      row_names_max_width = unit(10, "cm")  # Prevent label cutoff
    )
    
    return(ht)
  }, error = function(e) {
    cat(paste("Error creating heatmap:", e$message, "\n"))
    return(NULL)
  })
}

# =============================================================================
# SECTION 6: Compare Enrichment Results Between Datasets
# =============================================================================

cat("\n=============================================================================\n")
cat("COMPARING ENRICHMENT ACROSS DATASETS\n")
cat("=============================================================================\n\n")

# Function to compare pathways between datasets
compare_enrichment <- function(enrich1, enrich2, name1, name2, type = "KEGG") {
  
  if (is.null(enrich1[[type]]) || is.null(enrich2[[type]])) {
    cat(paste("Cannot compare", type, "- missing results\n"))
    return(NULL)
  }
  
  # Extract pathway IDs
  paths1 <- as.data.frame(enrich1[[type]])$ID
  paths2 <- as.data.frame(enrich2[[type]])$ID
  
  # Find common and unique pathways
  common <- intersect(paths1, paths2)
  unique1 <- setdiff(paths1, paths2)
  unique2 <- setdiff(paths2, paths1)
  
  cat(paste("\n", type, "Pathway Comparison:\n"))
  cat(paste("Common pathways:", length(common), "\n"))
  cat(paste("Unique to", name1, ":", length(unique1), "\n"))
  cat(paste("Unique to", name2, ":", length(unique2), "\n"))
  
  # Create comparison visualization
  comparison_data <- data.frame(
    Category = c("Common", paste("Unique to", name1), paste("Unique to", name2)),
    Count = c(length(common), length(unique1), length(unique2))
  )
  
  return(list(common = common, unique1 = unique1, unique2 = unique2, 
              data = comparison_data))
}

# Compare enrichment results if available
if (!is.null(enrichment_202220) && !is.null(enrichment_213323)) {
  kegg_comparison <- compare_enrichment(enrichment_202220, enrichment_213323,
                                       "GSE202220", "GSE213323", "KEGG")
  
  # Create Venn diagram-style comparison plot
  if (!is.null(kegg_comparison)) {
    pdf(file.path(plots_dir, "KEGG_pathway_comparison.pdf"), width = 8, height = 6)
    barplot(kegg_comparison$data$Count, 
            names.arg = kegg_comparison$data$Category,
            col = c("purple", "red", "blue"),
            main = "KEGG Pathway Comparison Between Datasets",
            ylab = "Number of Pathways",
            las = 2)
    dev.off()
  }
}

# =============================================================================
# SECTION 7: Session Info for Reproducibility
# =============================================================================

cat("\n=============================================================================\n")
cat("SAVING SESSION INFORMATION\n")
cat("=============================================================================\n\n")

# Capture session information (R version, package versions, etc.)
session_file <- file.path(output_dir, "session_info.txt")
sink(session_file)
cat("RNA-seq Analysis Session Information\n")
cat("=====================================\n\n")
cat(paste("Analysis Date:", Sys.time(), "\n\n"))
sessionInfo()
sink()

cat(paste("Session information saved to:", session_file, "\n"))

# =============================================================================
# SECTION 8: Generate Final Report Summary
# =============================================================================

cat("\n=============================================================================\n")
cat("GENERATING FINAL REPORT\n")
cat("=============================================================================\n\n")

# Create a comprehensive summary report
report_file <- file.path(output_dir, "Analysis_Report.txt")
sink(report_file)

cat("###############################################################################\n")
cat("# RNA-seq Differential Expression Analysis Report\n")
cat("# Generated:", as.character(Sys.time()), "\n")
cat("###############################################################################\n\n")

cat("DATASETS ANALYZED:\n")
cat("==================\n")
cat("1. GSE202220\n")
cat("2. GSE213323\n\n")

cat("ANALYSIS WORKFLOW:\n")
cat("==================\n")
cat("1. Raw count data download from GEO\n")
cat("2. Quality control and filtering\n")
cat("3. DESeq2 differential expression analysis\n")
cat("4. Multiple visualization (PCA, MA plots, volcano plots)\n")
cat("5. Pathway enrichment analysis (GO, KEGG, Reactome)\n")
cat("6. Gene Set Enrichment Analysis (GSEA)\n")
cat("7. Cross-dataset comparison\n\n")

cat("SUMMARY STATISTICS:\n")
cat("===================\n")
print(summary_stats)
cat("\n")

cat("OUTPUT FILES:\n")
cat("=============\n")
cat("Tables:\n")
cat("  - GSE202220_DEGs.txt: Full DESeq2 results for GSE202220\n")
cat("  - GSE213323_DEGs.txt: Full DESeq2 results for GSE213323\n")
cat("  - GSE202220_GSE213323_DEGAnalysis.xls: Combined results\n")
cat("  - *_enrichment.txt: Pathway enrichment results\n")
cat("  - *_GSEA_*.txt: Gene set enrichment results\n")
cat("  - DEG_Summary_Statistics.txt: Overview statistics\n\n")

cat("Plots:\n")
cat("  - *_PCA_plot.pdf: Principal component analysis\n")
cat("  - *_MA_plot.pdf: Log ratio vs abundance plots\n")
cat("  - *_volcano_plot.pdf: Volcano plots of DEGs\n")
cat("  - *_dotplot.pdf: Enrichment dotplots\n")
cat("  - *_GSEA_enrichment_plot.pdf: GSEA results\n")
cat("  - DEG_summary_barplot.pdf: Overall DEG summary\n\n")

cat("Quality Control:\n")
cat("  - *_QC_plots.pdf: Dispersion, library size, gene detection plots\n\n")

cat("INTERPRETATION GUIDELINES:\n")
cat("==========================\n")
cat("1. Significant DEGs: padj < 0.05\n")
cat("2. Fold change threshold: |log2FC| > 1 (2-fold change)\n")
cat("3. Blue points in volcano plots: downregulated genes\n")
cat("4. Red points in volcano plots: upregulated genes\n")
cat("5. Enrichment FDR < 0.05 considered significant\n\n")

cat("NEXT STEPS:\n")
cat("===========\n")
cat("1. Review volcano plots for genes of interest\n")
cat("2. Examine enriched pathways for biological insights\n")
cat("3. Validate key findings with qRT-PCR or other methods\n")
cat("4. Consider overlap analysis between datasets\n")
cat("5. Investigate top DEGs in literature\n\n")

cat("###############################################################################\n")

sink()

cat(paste("Final report saved to:", report_file, "\n"))

# =============================================================================
# SECTION 9: Create Example Pathway Heatmap (if enrichment results available)
# =============================================================================

cat("\n=============================================================================\n")
cat("CREATING EXAMPLE PATHWAY HEATMAPS\n")
cat("=============================================================================\n\n")

# Example: Create heatmap for top enriched pathway if KEGG results exist
if (!is.null(enrichment_202220$KEGG) && nrow(as.data.frame(enrichment_202220$KEGG)) > 0) {
  
  tryCatch({
    # Get top pathway
    top_pathway <- as.data.frame(enrichment_202220$KEGG)[1, ]
    
    cat(paste("Creating heatmap for top pathway:", top_pathway$Description, "\n"))
    
    # Create heatmap
    pdf(file.path(plots_dir, "Top_Pathway_Heatmap.pdf"), width = 10, height = 8)
    ht <- make_selected_pathway_heatmap(
      selectAnnotate = top_pathway,
      DEseqrbound = DEseqrbound,
      gene_col = "SYMBOL",
      comparison_col = "Subset_Comparison",
      value_col = "log2FoldChange",
      title = top_pathway$Description
    )
    if (!is.null(ht)) {
      draw(ht)
    }
    dev.off()
  }, error = function(e) {
    cat(paste("Could not create pathway heatmap:", e$message, "\n"))
  })
}

# =============================================================================
# SECTION 10: Identify Common DEGs Between Datasets
# =============================================================================

cat("\n=============================================================================\n")
cat("IDENTIFYING COMMON DEGs ACROSS DATASETS\n")
cat("=============================================================================\n\n")

# Find genes that are significant in both datasets
sig_202220 <- GSE202220_enrichment %>%
  filter(padj < 0.05, !is.na(SYMBOL)) %>%
  pull(SYMBOL)

sig_213323 <- GSE213323_enrichment %>%
  filter(padj < 0.05, !is.na(SYMBOL)) %>%
  pull(SYMBOL)

# Common DEGs
common_degs <- intersect(sig_202220, sig_213323)

cat(paste("Significant DEGs in GSE202220:", length(sig_202220), "\n"))
cat(paste("Significant DEGs in GSE213323:", length(sig_213323), "\n"))
cat(paste("Common DEGs across both datasets:", length(common_degs), "\n\n"))

# Save common DEGs
if (length(common_degs) > 0) {
  common_deg_data <- DEseqrbound %>%
    filter(SYMBOL %in% common_degs) %>%
    arrange(SYMBOL, Subset_Comparison)
  
  write.table(common_deg_data,
              file.path(tables_dir, "Common_DEGs_Both_Datasets.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Create heatmap of common DEGs
  if (length(common_degs) <= 100) {  # Only if manageable number
    pdf(file.path(plots_dir, "Common_DEGs_Heatmap.pdf"), width = 8, height = 12)
    
    common_wide <- DEseqrbound %>%
      filter(SYMBOL %in% common_degs) %>%
      select(SYMBOL, Subset_Comparison, log2FoldChange) %>%
      tidyr::pivot_wider(names_from = Subset_Comparison, 
                        values_from = log2FoldChange,
                        values_fn = mean)
    
    mat <- as.matrix(common_wide[, -1])
    rownames(mat) <- common_wide$SYMBOL
    
    ht <- Heatmap(
      mat,
      name = "log2FC",
      column_title = "Common DEGs Across Datasets",
      row_names_gp = gpar(fontsize = 6),
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
    )
    draw(ht)
    dev.off()
  }
}

# Create Venn diagram comparison
pdf(file.path(plots_dir, "DEG_Venn_Comparison.pdf"), width = 8, height = 6)
par(mar = c(2, 2, 2, 2))
# Create a simple representation
plot(1, type = "n", xlim = c(0, 10), ylim = c(0, 10), 
     xlab = "", ylab = "", axes = FALSE, main = "DEG Overlap Between Datasets")
# GSE202220 circle
symbols(3, 5, circles = 2, inches = FALSE, add = TRUE, fg = "red", lwd = 2)
text(2, 5, paste("GSE202220\n", length(sig_202220)), col = "red", cex = 1.2)
# GSE213323 circle
symbols(7, 5, circles = 2, inches = FALSE, add = TRUE, fg = "blue", lwd = 2)
text(8, 5, paste("GSE213323\n", length(sig_213323)), col = "blue", cex = 1.2)
# Overlap
text(5, 5, paste("Common\n", length(common_degs)), col = "purple", cex = 1.4, font = 2)
dev.off()

# =============================================================================
# FINAL SUMMARY
# =============================================================================

cat("\n###############################################################################\n")
cat("# ANALYSIS COMPLETE!\n")
cat("###############################################################################\n\n")

cat("OUTPUT DIRECTORY STRUCTURE:\n")
cat("===========================\n")
cat(paste(output_dir, "/\n"))
cat("  ├── tables/           # All result tables and enrichment data\n")
cat("  ├── plots/            # All visualization files\n")
cat("  ├── quality_control/  # QC plots and metrics\n")
cat("  ├── session_info.txt  # R session information\n")
cat("  └── Analysis_Report.txt # Comprehensive summary report\n\n")

cat("KEY FINDINGS:\n")
cat("=============\n")
cat(paste("- Total genes analyzed: ~", format(nrow(DEseqrbound)/2, big.mark = ","), "\n"))
cat(paste("- Significant DEGs (GSE202220):", summary_stats$Significant_DEGs[1], "\n"))
cat(paste("- Significant DEGs (GSE213323):", summary_stats$Significant_DEGs[2], "\n"))
cat(paste("- Common DEGs across both:", length(common_degs), "\n\n"))

