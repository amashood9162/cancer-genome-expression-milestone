# Load libraries
library(DESeq2)

# Load raw data
counts <- read.delim("data/TCGA-BRCA.star_counts.tsv.gz", row.names = 1)
clinical <- read.delim("data/TCGA-BRCA.clinical.tsv.gz")

# Reverse log2 transform to get raw counts
counts_raw <- round(2^counts - 1)

# Fix barcode formatting to match
clinical$sample_fixed <- gsub("-", ".", clinical$sample)

# Keep only Primary Tumor and Solid Tissue Normal samples
clinical_filtered <- clinical[clinical$sample_type.samples %in% c("Primary Tumor", "Solid Tissue Normal"), ]

# Keep only samples that have matching expression data
clinical_filtered <- clinical_filtered[clinical_filtered$sample_fixed %in% colnames(counts_raw), ]

# Match counts columns to clinical row order
counts_filtered <- counts_raw[, clinical_filtered$sample_fixed]

# Verify alignment
stopifnot(all(colnames(counts_filtered) == clinical_filtered$sample_fixed))

# Build DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = clinical_filtered,
  design = ~ sample_type.samples
)

dds$sample_type.samples <- relevel(dds$sample_type.samples, ref = "Solid Tissue Normal")

dds
saveRDS(dds, "data/dds_raw.rds")