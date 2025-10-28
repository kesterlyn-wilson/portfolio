# FASTQ QC Summary

**Description:**  
This lightweight Python script parses FASTQ quality metric files and outputs a summary table of key run metrics, including read counts, GC content, average quality scores, and sequence length distribution. It’s designed for quick QC checks of raw sequencing data prior to downstream analysis.

**Use Case**  
The script can be used by NGS support teams, bioinformaticians, or lab scientists who want to:
- Validate sequencing run quality
- Compare read metrics across samples
- Identify low-quality or biased libraries before alignment

---

### 🚀 Features
- Parses standard FASTQ or text-based QC outputs  
- Computes key metrics (mean Q-score, read length, GC%)  
- Exports results as a clean CSV summary  
- Compatible with outputs from tools like FastQC, SeqKit, or custom log files

---

### 🧰 Dependencies
pandas
biopython

### Install them:
pip install -r requirements.txt