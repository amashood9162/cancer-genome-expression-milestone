library(DESeq2)
library(ggplot2)

# Load the analyzed dataset
dds <- readRDS("data/dds_analyzed.rds")
res_df <- read.csv("results/full_DE_results.csv")

# ---- PLOT 1: PCA Plot ----
# Variance-stabilizing transform: makes data suitable for PCA
# (raw counts have unequal variance across expression levels, VST corrects this)
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "sample_type.samples", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p1 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = sample_type.samples)) +
  geom_point(size = 2, alpha = 0.7) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA: Tumor vs Normal Samples") +
  theme_minimal() +
  scale_color_manual(values = c("Primary Tumor" = "firebrick", "Solid Tissue Normal" = "steelblue"))

ggsave("results/pca_plot.png", p1, width = 7, height = 5, dpi = 300)

# ---- PLOT 2: Volcano Plot ----
res_df$significant <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                             "Significant", "Not Significant")

p2 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c("Significant" = "firebrick", "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  xlab("Log2 Fold Change (Tumor vs Normal)") +
  ylab("-Log10 Adjusted P-value") +
  ggtitle("Volcano Plot: Differential Gene Expression") +
  theme_minimal()

ggsave("results/volcano_plot.png", p2, width = 7, height = 6, dpi = 300)

# ---- PLOT 3: Heatmap of top 30 significant genes ----
library(pheatmap)

sig_genes <- read.csv("results/significant_genes.csv")
top30 <- head(sig_genes[order(sig_genes$padj), ], 30)

vsd_mat <- assay(vsd)[top30$gene_id, ]

# Annotation showing which samples are tumor vs normal
annotation_col <- data.frame(SampleType = colData(dds)$sample_type.samples)
rownames(annotation_col) <- colnames(vsd_mat)

png("results/heatmap_top30.png", width = 1200, height = 1000, res = 150)
pheatmap(vsd_mat,
         scale = "row",
         show_colnames = FALSE,
         annotation_col = annotation_col,
         main = "Top 30 Significant Genes: Tumor vs Normal")
dev.off()

cat("All plots saved to results/ folder\n")