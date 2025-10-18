#!/usr/bin/env python3
"""
FASTQ QC Summary Script
Author: Kesterlyn Wilson
Date: 2025-10
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO

def calculate_gc(seq):
    """Calculate GC content for a given sequence."""
    g = seq.count('G')
    c = seq.count('C')
    gc_content = (g + c) / len(seq) * 100 if len(seq) > 0 else 0
    return gc_content

def summarize_fastq(file_path):
    """Compute basic QC metrics for a FASTQ file."""
    total_reads = 0
    total_bases = 0
    total_gc = 0
    total_q = 0
    seq_lengths = []

    for record in SeqIO.parse(file_path, "fastq"):
        seq = str(record.seq)
        q_scores = record.letter_annotations["phred_quality"]

        total_reads += 1
        total_bases += len(seq)
        seq_lengths.append(len(seq))
        total_gc += calculate_gc(seq)
        total_q += sum(q_scores) / len(q_scores)

    avg_qscore = total_q / total_reads if total_reads else 0
    avg_gc = total_gc / total_reads if total_reads else 0
    avg_len = total_bases / total_reads if total_reads else 0

    sample_id = os.path.basename(file_path).split(".")[0]

    return {
        "Sample_ID": sample_id,
        "Total_Reads": total_reads,
        "Avg_QScore": round(avg_qscore, 2),
        "GC_Content": round(avg_gc, 2),
        "Mean_Read_Length": round(avg_len, 1),
    }

def plot_qc_summary(df):
    """Generate basic QC summary plots."""
    # Plot 1: Mean read length per sample
    plt.figure(figsize=(6,4))
    plt.bar(df['Sample_ID'], df['Mean_Read_Length'])
    plt.title("Mean Read Length per Sample")
    plt.xlabel("Sample ID")
    plt.ylabel("Mean Read Length (bp)")
    plt.tight_layout()
    plt.savefig("output/mean_read_length.png")
    plt.close()

    # Plot 2: GC content per sample
    plt.figure(figsize=(6,4))
    plt.bar(df['Sample_ID'], df['GC_Content'])
    plt.title("Average GC Content per Sample")
    plt.xlabel("Sample ID")
    plt.ylabel("GC %")
    plt.tight_layout()
    plt.savefig("output/gc_content.png")
    plt.close()

def main(input_dir):
    """Main function to process all FASTQ files and generate QC reports."""
    summaries = []
    for file in os.listdir(input_dir):
        if file.endswith(".fastq") or file.endswith(".fq"):
            path = os.path.join(input_dir, file)
            print(f"Processing {file}...")
            summary = summarize_fastq(path)
            summaries.append(summary)

    if not summaries:
        print("No FASTQ files found in the input directory.")
        return

    df = pd.DataFrame(summaries)
    os.makedirs("output", exist_ok=True)
    df.to_csv("output/fastq_qc_summary.csv", index=False)
    plot_qc_summary(df)
    print("✅ QC summary and plots saved to output/")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python fastq_qc_summary.py <input_directory>")
        sys.exit(1)
    
    input_directory = sys.argv[1]
    
    if not os.path.isdir(input_directory):
        print(f"Error: {input_directory} is not a valid directory")
        sys.exit(1)
    
    main(input_directory)