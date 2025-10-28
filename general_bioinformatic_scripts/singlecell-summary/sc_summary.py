#!/usr/bin/env python3
"""
Single-Cell Expression Summary
Author: Kesterlyn Wilson
Date: 2025-10
"""
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def summarize_singlecell(file_path):
    """Analyze single-cell expression matrix and generate QC metrics."""
    print(f"Processing {os.path.basename(file_path)}...")
    
    try:
        # Read expression matrix (genes x cells)
        df = pd.read_csv(file_path, index_col=0)
        print(f"  Matrix shape: {df.shape[0]} genes × {df.shape[1]} cells")
        
    except Exception as e:
        print(f"  ❌ Error reading file: {e}")
        return
    
    # Basic stats per cell
    counts_per_cell = df.sum(axis=0)
    genes_per_cell = (df > 0).sum(axis=0)
    
    # Basic stats per gene
    mean_expression = df.mean(axis=1)
    cells_per_gene = (df > 0).sum(axis=1)
    
    # Create output directory
    os.makedirs("output", exist_ok=True)
    
    # Save per-cell summary
    cell_summary = pd.DataFrame({
        "Total_Counts": counts_per_cell,
        "Genes_Detected": genes_per_cell
    })
    cell_summary.to_csv("output/sc_summary_per_cell.csv")
    
    # Save per-gene summary
    gene_summary = pd.DataFrame({
        "Mean_Expression": mean_expression,
        "Cells_Expressing": cells_per_gene,
        "Detection_Rate": (cells_per_gene / df.shape[1] * 100).round(2)
    })
    gene_summary.to_csv("output/sc_summary_per_gene.csv")
    
    # Generate plots
    generate_qc_plots(counts_per_cell, genes_per_cell)
    
    # Print summary statistics
    print(f"  Summary Statistics:")
    print(f"    Mean counts/cell: {counts_per_cell.mean():.0f}")
    print(f"    Median counts/cell: {counts_per_cell.median():.0f}")
    print(f"    Mean genes/cell: {genes_per_cell.mean():.0f}")
    print(f"    Median genes/cell: {genes_per_cell.median():.0f}")
    print("  ✅ Analysis complete!")

def generate_qc_plots(counts_per_cell, genes_per_cell):
    """Generate QC visualization plots."""
    
    # Determine optimal bin count (Sturges' rule)
    n_bins_counts = min(50, int(np.ceil(np.log2(len(counts_per_cell)) + 1)))
    n_bins_genes = min(50, int(np.ceil(np.log2(len(genes_per_cell)) + 1)))
    
    # Plot 1: Total counts per cell histogram
    plt.figure(figsize=(8, 5))
    plt.hist(counts_per_cell, bins=n_bins_counts, edgecolor='black', alpha=0.7)
    plt.axvline(counts_per_cell.median(), color='red', linestyle='--', 
                linewidth=2, label=f'Median: {counts_per_cell.median():.0f}')
    plt.title("Total UMI Counts per Cell", fontsize=14, fontweight='bold')
    plt.xlabel("UMI Counts", fontsize=12)
    plt.ylabel("Number of Cells", fontsize=12)
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig("output/total_counts_hist.png", dpi=150)
    plt.close()
    
    # Plot 2: Genes detected per cell histogram
    plt.figure(figsize=(8, 5))
    plt.hist(genes_per_cell, bins=n_bins_genes, edgecolor='black', alpha=0.7)
    plt.axvline(genes_per_cell.median(), color='red', linestyle='--', 
                linewidth=2, label=f'Median: {genes_per_cell.median():.0f}')
    plt.title("Genes Detected per Cell", fontsize=14, fontweight='bold')
    plt.xlabel("Number of Genes", fontsize=12)
    plt.ylabel("Number of Cells", fontsize=12)
    plt.legend()
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig("output/genes_per_cell_hist.png", dpi=150)
    plt.close()
    
    # Plot 3: Scatter plot - counts vs genes detected
    plt.figure(figsize=(8, 6))
    plt.scatter(counts_per_cell, genes_per_cell, alpha=0.5, s=10)
    plt.title("UMI Counts vs Genes Detected", fontsize=14, fontweight='bold')
    plt.xlabel("Total UMI Counts", fontsize=12)
    plt.ylabel("Genes Detected", fontsize=12)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("output/counts_vs_genes_scatter.png", dpi=150)
    plt.close()

def main(input_dir):
    """Process all CSV files in the input directory."""
    
    # Validate input directory
    if not os.path.isdir(input_dir):
        print(f"❌ Error: '{input_dir}' is not a valid directory")
        sys.exit(1)
    
    # Find CSV files
    csv_files = [f for f in os.listdir(input_dir) if f.endswith(".csv")]
    
    if not csv_files:
        print(f"❌ No CSV files found in '{input_dir}'")
        sys.exit(1)
    
    print(f"Found {len(csv_files)} CSV file(s)\n")
    
    # Process each file
    for file in csv_files:
        summarize_singlecell(os.path.join(input_dir, file))
        print()
    
    print("=" * 60)
    print("📊 All outputs saved to output/")
    print("  - sc_summary_per_cell.csv")
    print("  - sc_summary_per_gene.csv")
    print("  - total_counts_hist.png")
    print("  - genes_per_cell_hist.png")
    print("  - counts_vs_genes_scatter.png")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python singlecell_summary.py <input_directory>")
        print("\nExpected input: CSV file with genes as rows, cells as columns")
        sys.exit(1)
    
    main(sys.argv[1])