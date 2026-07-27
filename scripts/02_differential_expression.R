library(DESeq2)

# Load the prepared dataset (in case this is a fresh session)
dds <- readRDS("data/dds_raw.rds")

# Run the actual differential expression test
# This fits DESeq2's statistical model across all ~60,000 genes
# With 1,219 samples this can take 10-20+ minutes - be patient
dds <- DESeq(dds)

# Extract results: log2 fold change, p-value, adjusted p-value per gene
# alpha = 0.05 tells DESeq2 our significance threshold for summary purposes
res <- results(dds, alpha = 0.05)

# Sort genes by adjusted p-value (most statistically confident first)
res_ordered <- res[order(res$padj), ]

# Quick summary: how many genes are significant, up vs down
summary(res)

# Convert to a data frame and add gene IDs as a real column
res_df <- as.data.frame(res_ordered)
res_df$gene_id <- rownames(res_df)

# Filter to significant, meaningfully-changed genes:
# adjusted p-value < 0.05 AND at least a 2-fold change (log2FC of 1 = 2x)
sig_genes <- res_df[!is.na(res_df$padj) & 
                      res_df$padj < 0.05 & 
                      abs(res_df$log2FoldChange) > 1, ]

# How many significant genes did we find?
nrow(sig_genes)

# Save full results and significant-only results as CSV files
write.csv(res_df, "results/full_DE_results.csv", row.names = FALSE)
write.csv(sig_genes, "results/significant_genes.csv", row.names = FALSE)

# Save the dds object with results baked in, for later use (e.g. plotting)
saveRDS(dds, "data/dds_analyzed.rds")