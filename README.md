# neutrophil-rnaseq-demo
# Neutrophil RNA-seq Demo

A tiny project to practice basic bulk RNA-seq analysis on an example count matrix.  
The goal is to simulate differential gene expression between control and lesional skin (or neutrophils vs control) and visualize inflammatory signatures.

### Biological Context

This toy dataset mimics the transcriptomic profiling of inflammatory lesions:
- cytokines and chemokines (IL1A, IL6, CXCL8, TNF, TSLP, CCL5, etc.)
- metalloproteases (MMP2, MMP8, MMP9)
- neutrophil mediators (MPO, ELANE, LTF, LCN2, CAMP, etc.)

We expect inflammatory mediators and neutrophil activation markers to be upregulated in lesional tissue.

### What this mini-project does

- Loads a count matrix (gene × sample)
- Defines sample conditions (Control vs Lesion)
- Computes log2 fold changes
- Highlights inflammatory genes
- Plots:
  - volcano plot
  - small heatmap of selected genes

### Files

- `counts_matrix.csv` — toy gene counts
- `sample_metadata.csv` — sample → condition mapping
- `analyze_rnaseq.py` — code for differential expression and plotting

### How to run

Install dependencies:

```bash
pip install pandas numpy matplotlib seaborn
