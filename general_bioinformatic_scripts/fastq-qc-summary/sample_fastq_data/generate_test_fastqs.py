#!/usr/bin/env python3
"""
Generate random test FASTQ files for QC testing
"""

import random
import os

def generate_sequence(length):
    """Generate a random DNA sequence."""
    return ''.join(random.choices('ACGT', k=length))

def generate_quality_scores(length, avg_quality=30):
    """Generate Phred quality scores as ASCII characters."""
    # Phred+33 encoding: Q score + 33 = ASCII character
    # Q30 = '?' (ASCII 63), Q40 = 'I' (ASCII 73)
    qualities = []
    for _ in range(length):
        q = max(0, min(40, int(random.gauss(avg_quality, 5))))
        qualities.append(chr(q + 33))
    return ''.join(qualities)

def write_fastq_file(filename, num_reads=1000, read_length=150, avg_quality=30):
    """Write a FASTQ file with specified parameters."""
    with open(filename, 'w') as f:
        for i in range(1, num_reads + 1):
            # Vary read length slightly
            length = read_length + random.randint(-10, 10)
            seq = generate_sequence(length)
            qual = generate_quality_scores(length, avg_quality)
            
            # Write FASTQ record (4 lines per read)
            f.write(f"@READ_{i}\n")
            f.write(f"{seq}\n")
            f.write("+\n")
            f.write(f"{qual}\n")

def main():
    """Generate multiple test FASTQ files with different characteristics."""
    
    # Create test_data directory
    os.makedirs("test_data", exist_ok=True)
    
    # Sample 1: High quality, normal GC content
    print("Generating sample1.fastq...")
    write_fastq_file("test_data/sample1.fastq", 
                     num_reads=500, 
                     read_length=150, 
                     avg_quality=35)
    
    # Sample 2: Medium quality, slightly longer reads
    print("Generating sample2.fastq...")
    write_fastq_file("test_data/sample2.fastq", 
                     num_reads=500, 
                     read_length=200, 
                     avg_quality=28)
    
    # Sample 3: Lower quality, shorter reads
    print("Generating sample3.fastq...")
    write_fastq_file("test_data/sample3.fastq", 
                     num_reads=500, 
                     read_length=100, 
                     avg_quality=22)
    
    print("\n✅ Test FASTQ files created in test_data/")
    print("   - sample1.fastq (500 reads, ~150bp, high quality)")
    print("   - sample2.fastq (500 reads, ~200bp, medium quality)")
    print("   - sample3.fastq (500 reads, ~100bp, lower quality)")
    print("\nRun QC script with: python fastq_qc_summary.py test_data/")

if __name__ == "__main__":
    main()