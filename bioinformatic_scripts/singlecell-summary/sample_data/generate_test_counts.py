#!/usr/bin/env python3
"""
Generate test single-cell expression count matrix
Simulates realistic UMI count data
"""

import numpy as np
import pandas as pd
import os

def generate_count_matrix(n_genes=2000, n_cells=500, seed=None):
    """
    Generate a realistic single-cell count matrix.
    
    Parameters:
    -----------
    n_genes : int
        Number of genes (rows)
    n_cells : int
        Number of cells (columns)
    seed : int, optional
        Random seed for reproducibility
    
    Returns:
    --------
    pd.DataFrame
        Count matrix with genes as rows, cells as columns
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Initialize empty matrix
    counts = np.zeros((n_genes, n_cells), dtype=int)
    
    # Simulate different cell types with varying quality
    for i in range(n_cells):
        # Cell quality varies (some cells have higher counts/more genes)
        if i < n_cells * 0.7:  # 70% high-quality cells
            mean_counts = np.random.uniform(2000, 5000)
            dropout_rate = 0.85  # 85% of genes are zeros
        elif i < n_cells * 0.9:  # 20% medium-quality cells
            mean_counts = np.random.uniform(1000, 2000)
            dropout_rate = 0.90
        else:  # 10% low-quality cells
            mean_counts = np.random.uniform(500, 1000)
            dropout_rate = 0.95
        
        # Generate expression for this cell
        for j in range(n_genes):
            # Dropout: most genes are not expressed
            if np.random.random() > dropout_rate:
                # Use negative binomial for realistic count distribution
                # (overdispersed compared to Poisson)
                mean_expr = np.random.exponential(mean_counts / (1 - dropout_rate) / 100)
                count = np.random.poisson(mean_expr)
                counts[j, i] = count
    
    # Create gene and cell names
    gene_names = [f"Gene_{i+1}" for i in range(n_genes)]
    cell_names = [f"Cell_{i+1}" for i in range(n_cells)]
    
    # Create DataFrame
    df = pd.DataFrame(counts, index=gene_names, columns=cell_names)
    
    return df

def generate_marker_genes(df, n_markers=50):
    """
    Add some highly expressed marker genes to make data more realistic.
    """
    df_copy = df.copy()
    
    # Select random genes to be markers
    marker_indices = np.random.choice(df.shape[0], size=n_markers, replace=False)
    
    for idx in marker_indices:
        # Markers are expressed in most cells at higher levels
        for col in df_copy.columns:
            if np.random.random() > 0.3:  # 70% of cells express marker
                df_copy.iloc[idx][col] = np.random.poisson(50) + 10
    
    return df_copy

def print_matrix_stats(df, name="Matrix"):
    """Print summary statistics of the matrix."""
    print(f"\n{name} Statistics:")
    print(f"  Shape: {df.shape[0]} genes × {df.shape[1]} cells")
    print(f"  Total UMIs: {df.sum().sum():,}")
    print(f"  Mean counts/cell: {df.sum(axis=0).mean():.0f}")
    print(f"  Median counts/cell: {df.sum(axis=0).median():.0f}")
    print(f"  Mean genes/cell: {(df > 0).sum(axis=0).mean():.0f}")
    print(f"  Median genes/cell: {(df > 0).sum(axis=0).median():.0f}")
    print(f"  Sparsity: {(df == 0).sum().sum() / df.size * 100:.1f}%")

def main():
    """Generate test single-cell count matrices."""
    
    # Create test_data directory
    os.makedirs("test_sc_data", exist_ok=True)
    
    print("Generating test single-cell count matrices...\n")
    print("=" * 60)
    
    # Matrix 1: Standard dataset
    print("\n📊 Generating dataset1.csv (standard quality)...")
    df1 = generate_count_matrix(n_genes=2000, n_cells=500, seed=42)
    df1 = generate_marker_genes(df1, n_markers=50)
    df1.to_csv("test_sc_data/dataset1.csv")
    print_matrix_stats(df1, "Dataset 1")
    
    # Matrix 2: Smaller dataset with lower quality
    print("\n📊 Generating dataset2.csv (lower quality)...")
    df2 = generate_count_matrix(n_genes=1500, n_cells=300, seed=43)
    # Don't add as many markers for this one
    df2 = generate_marker_genes(df2, n_markers=20)
    df2.to_csv("test_sc_data/dataset2.csv")
    print_matrix_stats(df2, "Dataset 2")
    
    # Matrix 3: High-quality dataset
    print("\n📊 Generating dataset3.csv (high quality)...")
    np.random.seed(44)
    df3 = generate_count_matrix(n_genes=2500, n_cells=400, seed=44)
    df3 = generate_marker_genes(df3, n_markers=80)
    df3.to_csv("test_sc_data/dataset3.csv")
    print_matrix_stats(df3, "Dataset 3")
    
    print("\n" + "=" * 60)
    print("✅ Test count matrices created in test_sc_data/")
    print("\nFiles created:")
    print("  - dataset1.csv (2000 genes × 500 cells)")
    print("  - dataset2.csv (1500 genes × 300 cells)")
    print("  - dataset3.csv (2500 genes × 400 cells)")
    print("\nRun analysis with:")
    print("  python singlecell_summary.py test_sc_data/")

if __name__ == "__main__":
    main()