# FASTQ QC Summary

**Description:**  
Very lightweight Python script to parse FASTQ quality metric files w/ a summary table of key run metrics (read counts, GC content, average quality scores, and sequence length distribution). Conduct quick QC checks of raw sequencing data prior to downstream analysis.
- Parses standard FASTQ or text-based QC outputs  
- Computes key quality metrics (mean Q-score, read length, GC%)  
- Exports results in CSV summary  
- Compatible with outputs from tools like FastQC, SeqKit, or custom log files

---

### Dependencies
pandas
biopython

### Install them:
pip install -r requirements.txt
