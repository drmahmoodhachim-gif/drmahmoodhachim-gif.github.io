# =============================================================================
# RNA-seq Tutorial - EVERYTHING IN R (No Python needed!)
# =============================================================================
# For PhD students with no programming background.
# 
# HOW TO RUN:
#   Option 1: Double-click "RUN_ANALYSIS.bat" (Windows)
#   Option 2: Open RStudio, open this file, click "Source" button
#   Option 3: In R console: source("run_all_in_R.R")
#
# First time: Install packages (run once in R):
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("DESeq2","ggplot2","pheatmap","RColorBrewer"))
# =============================================================================

cat("\n")
cat("========================================\n")
cat("  RNA-seq Pipeline - All in R\n")
cat("  (No Python required!)\n")
cat("========================================\n\n")

# --- Step 1: Generate synthetic count matrix (simulates RNA-seq data) ---
cat("[Step 1] Creating synthetic RNA-seq data...\n")

set.seed(42)
genes <- paste0("GENE_", sprintf("%03d", 0:49))
samples <- c(paste0("CTRL_", sprintf("%02d", 1:4)), paste0("TREAT_", sprintf("%02d", 1:4)))

# Known DEGs: GENE_000-004 (up in treatment), GENE_045-049 (down in treatment)
DEG_UP   <- 1:5      # gene indices 1-5
DEG_DOWN <- 46:50    # gene indices 46-50

base <- 100
count_matrix <- matrix(0, nrow = 50, ncol = 8)

for (i in 1:50) {
  if (i %in% DEG_UP) {
    ctrl_vals <- base + sample(0:20, 4, replace = TRUE)
    treat_vals <- round((base + sample(0:20, 4, replace = TRUE)) * 3)
  } else if (i %in% DEG_DOWN) {
    ctrl_vals <- round((base + sample(0:20, 4, replace = TRUE)) * 3)
    treat_vals <- base + sample(0:20, 4, replace = TRUE)
  } else {
    ctrl_vals <- pmax(1, base + sample(-10:30, 4, replace = TRUE))
    treat_vals <- pmax(1, base + sample(-10:30, 4, replace = TRUE))
  }
  count_matrix[i, ] <- c(ctrl_vals, treat_vals)
}

count_matrix <- count_matrix + sample(-3:3, 400, replace = TRUE)
count_matrix <- pmax(1, count_matrix)

rownames(count_matrix) <- genes
colnames(count_matrix) <- samples

dir.create("04_counts", showWarnings = FALSE)
write.csv(count_matrix, "04_counts/count_matrix.csv")
cat("   Created: 04_counts/count_matrix.csv\n\n")

# --- Step 2: DEG Analysis ---
cat("[Step 2] Differential expression analysis (DESeq2)...\n")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

counts <- round(count_matrix)
metadata <- data.frame(
  sample = colnames(counts),
  condition = c(rep("Control", 4), rep("Treatment", 4)),
  row.names = colnames(counts)
)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Treatment", "Control"))

res_ordered <- res[order(res$padj), ]
res_df <- as.data.frame(res_ordered)
res_df$gene <- rownames(res_df)

sig <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Yes", "No")
res_df$significant[is.na(res_df$padj)] <- "No"

dir.create("05_deg_results", showWarnings = FALSE)
write.csv(res_df, "05_deg_results/all_genes.csv", row.names = FALSE)
write.csv(sig, "05_deg_results/significant_deg.csv", row.names = FALSE)

cat("   Significant DEGs:", nrow(sig), "\n\n")

# --- Step 3: Volcano plot ---
cat("[Step 3] Creating visualizations...\n")

p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj + 1e-100))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("No" = "gray70", "Yes" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  theme_minimal() +
  labs(title = "Volcano Plot - Treatment vs Control",
       x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
  theme(legend.position = "top")
ggsave("05_deg_results/volcano_plot.pdf", p_volcano, width = 6, height = 5)
ggsave("05_deg_results/volcano_plot.png", p_volcano, width = 6, height = 5, dpi = 150)

# --- Heatmap ---
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
if (nrow(sig) > 0) {
  top_genes <- rownames(sig)[1:min(20, nrow(sig))]
  mat_heat <- mat[top_genes, ]
  mat_scaled <- t(scale(t(mat_heat)))
  pdf("05_deg_results/heatmap_top_deg.pdf", width = 6, height = 6)
  pheatmap(mat_scaled, scale = "none",
           color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
           main = "Top DEGs - Z-scored expression", border_color = NA, fontsize_row = 8)
  dev.off()
}

# --- MA plot ---
p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("No" = "gray70", "Yes" = "red")) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  theme_minimal() +
  labs(title = "MA Plot - Treatment vs Control",
       x = "log10(mean expression)", y = "log2 Fold Change")
ggsave("05_deg_results/ma_plot.pdf", p_ma, width = 6, height = 5)

# --- Step 4: Systems biology ---
cat("[Step 4] Systems biology analysis...\n")

dir.create("06_systems_biology", showWarnings = FALSE)

if (nrow(sig) > 0) {
  deg_genes <- sig$gene
  fake_pathways <- list(
    "Cell cycle" = paste0("GENE_", sprintf("%03d", 0:9)),
    "Apoptosis" = paste0("GENE_", sprintf("%03d", 10:19)),
    "Metabolism" = paste0("GENE_", sprintf("%03d", 20:29)),
    "Signaling" = paste0("GENE_", sprintf("%03d", 30:39)),
    "Stress response" = paste0("GENE_", sprintf("%03d", 40:49))
  )
  pathway_enrichment <- data.frame(
    pathway = names(fake_pathways),
    deg_count = sapply(fake_pathways, function(p) sum(deg_genes %in% p)),
    total_genes = sapply(fake_pathways, length),
    stringsAsFactors = FALSE
  )
  pathway_enrichment$pct <- 100 * pathway_enrichment$deg_count / pathway_enrichment$total_genes

  p1 <- ggplot(pathway_enrichment, aes(x = reorder(pathway, deg_count), y = deg_count, fill = pathway)) +
    geom_col() + coord_flip() + theme_minimal() +
    labs(title = "DEGs by Biological Pathway", x = "", y = "Number of DEGs") +
    theme(legend.position = "none")
  ggsave("06_systems_biology/pathway_barplot.pdf", p1, width = 6, height = 4)
  ggsave("06_systems_biology/pathway_barplot.png", p1, width = 6, height = 4, dpi = 150)

  sig$regulation <- ifelse(sig$log2FoldChange > 0, "Upregulated", "Downregulated")
  reg_summary <- as.data.frame(table(sig$regulation))
  p2 <- ggplot(reg_summary, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_col() + theme_minimal() +
    labs(title = "DEG Regulation Summary", x = "", y = "Count") +
    scale_fill_brewer(palette = "Set2") + theme(legend.position = "none")
  ggsave("06_systems_biology/regulation_summary.pdf", p2, width = 5, height = 4)

  write.csv(pathway_enrichment, "06_systems_biology/pathway_enrichment.csv", row.names = FALSE)
  write.csv(reg_summary, "06_systems_biology/regulation_summary.csv", row.names = FALSE)
}

# --- Done! ---
cat("\n")
cat("========================================\n")
cat("  DONE! All results saved.\n")
cat("========================================\n")
cat("\nYour results:\n")
cat("  - 05_deg_results/volcano_plot.png     (open this to see your volcano plot!)\n")
cat("  - 05_deg_results/ma_plot.pdf\n")
cat("  - 05_deg_results/heatmap_top_deg.pdf\n")
cat("  - 05_deg_results/significant_deg.csv   (list of significant genes)\n")
cat("  - 06_systems_biology/pathway_barplot.png\n")
cat("\n")
cat("Open the PNG files by double-clicking them.\n")
cat("\n")
