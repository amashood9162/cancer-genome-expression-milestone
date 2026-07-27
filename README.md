# Cancer Gene Expression Analysis: TCGA Breast Cancer (BRCA)

## Overview
Differential gene expression analysis comparing tumor vs. normal breast tissue 
using TCGA-BRCA RNA-seq data (1,219 samples: 1,106 tumor, 113 normal), 
identifying candidate genes with therapeutic relevance.

## Methods
- **Data source:** UCSC Xena (GDC TCGA Breast Cancer cohort), STAR-Counts RNA-seq
- **Analysis:** DESeq2 (R/Bioconductor) for differential expression
- **Significance criteria:** adjusted p-value < 0.05, |log2FC| > 1

## Key Results
- 28,943 genes significantly differentially expressed (18,299 up, 10,644 down)
- 12,173 genes pass strict significance + fold-change filtering
- PCA confirms clear separation between tumor and normal samples at the 
  transcriptome level
- Top significant genes visualized via volcano plot and heatmap
-PCA analysis confirmed strong separation between tumor and normal transcriptomes,
validating the pipeline's biological accuracy before proceeding to differential expression testing.

## Files
- `scripts/01_data_prep.R` — data loading, cleaning, DESeq2 dataset construction
- `scripts/02_differential_expression.R` — DESeq2 statistical testing
- `scripts/03_visualization.R` — PCA, volcano plot, heatmap generation
- `results/` — output CSVs and plots

## Limitations & Next Steps
- Model does not yet control for covariates (e.g., tumor stage, age)
- Next: pathway enrichment analysis and cross-referencing candidate genes 
  against drug-gene interaction databases (e.g., DGIdb) to prioritize 
  plausible therapeutic targets from the 12,173-gene candidate list

## Tools
R, DESeq2, ggplot2, pheatmap