# RNA-seq Complete Manual: From Design to Analysis

**Step-by-step guide with R code you can copy and paste**

*For PhD students – designed like a ChatGPT tutorial: explanation → code → run*

---

## How to Use This Manual

1. Read each section in order.
2. Copy the R code blocks and paste them into R or RStudio.
3. Run the code and check the output.
4. Use the explanations to understand what each part does.

---

# Part 1: Designing Your RNA-seq Experiment

## Page 1.1 – Why Design Matters

Before any sequencing or analysis, you need a clear experimental design. Poor design leads to results you cannot interpret or publish.

**Key questions:**
- What are you comparing? (e.g., healthy vs diseased, control vs treatment)
- How many biological replicates? (at least 3 per condition, preferably more)
- How many conditions/groups?

---

## Page 1.2 – Basic Design Principles

| Principle | What it means | Example |
|-----------|----------------|---------|
| **Biological replicates** | Different animals, patients, or cell passages – not just technical repeats | 4 mice in control, 4 in treatment |
| **Randomization** | Randomly assign samples to groups to avoid bias | Don't run all controls on Monday and all treatments on Tuesday |
| **Balance** | Same number of replicates per condition | 4 vs 4, not 3 vs 6 |
| **Blocking** | If you know a source of variation (batch, sex), account for it in the design | Include "batch" or "sex" in your design formula |

---

## Page 1.3 – Example Designs

**Simple two-group comparison (most common):**
```
Control:    Sample 1, 2, 3, 4
Treatment:  Sample 5, 6, 7, 8
```

**With a blocking factor (e.g., batch):**
```
Batch 1: Control A, Control B, Treatment A, Treatment B
Batch 2: Control C, Control D, Treatment C, Treatment D
```

**Our tutorial uses:** 4 Control + 4 Treatment (simple two-group design).

---

# Part 2: Understanding the Data – Count Matrix

## Page 2.1 – What Is a Count Matrix?

After RNA-seq, you get a **count matrix**: a table where
- **Rows** = genes
- **Columns** = samples
- **Values** = number of reads mapping to each gene in each sample

**Example:**
```
         CTRL_01  CTRL_02  TREAT_01  TREAT_02
GENE_A      120      115       350       340
GENE_B       85       90        30        28
```

GENE_A is higher in treatment; GENE_B is lower. That’s what we want to find statistically.

---

## Page 2.2 – Create a Count Matrix in R (Our Tutorial Data)

Run this in R to build a small synthetic count matrix. Copy and paste the whole block.

```r
# ============================================
# Step: Create synthetic count matrix
# ============================================
# We simulate RNA-seq counts for 50 genes,
# 4 Control + 4 Treatment samples.
# Some genes are designed to be different (DEGs).
# ============================================

set.seed(42)
genes <- paste0("GENE_", sprintf("%03d", 0:49))
samples <- c(paste0("CTRL_", sprintf("%02d", 1:4)), 
             paste0("TREAT_", sprintf("%02d", 1:4)))

# DEGs: indices 1-5 UP in treatment, 46-50 DOWN
DEG_UP   <- 1:5
DEG_DOWN <- 46:50

base <- 100
count_matrix <- matrix(0, nrow = 50, ncol = 8)

for (i in 1:50) {
  if (i %in% DEG_UP) {
    ctrl_vals  <- base + sample(0:20, 4, replace = TRUE)
    treat_vals <- round((base + sample(0:20, 4, replace = TRUE)) * 3)
  } else if (i %in% DEG_DOWN) {
    ctrl_vals  <- round((base + sample(0:20, 4, replace = TRUE)) * 3)
    treat_vals <- base + sample(0:20, 4, replace = TRUE)
  } else {
    ctrl_vals  <- pmax(1, base + sample(-10:30, 4, replace = TRUE))
    treat_vals <- pmax(1, base + sample(-10:30, 4, replace = TRUE))
  }
  count_matrix[i, ] <- c(ctrl_vals, treat_vals)
}

count_matrix <- pmax(1, count_matrix + sample(-3:3, 400, replace = TRUE))
rownames(count_matrix) <- genes
colnames(count_matrix) <- samples

dir.create("04_counts", showWarnings = FALSE)
write.csv(count_matrix, "04_counts/count_matrix.csv")
cat("Saved: 04_counts/count_matrix.csv\n")
```

**What this does:**
- `set.seed(42)` – makes results reproducible
- `genes` – 50 gene names (GENE_000 to GENE_049)
- `samples` – 8 sample names (CTRL_01–04, TREAT_01–04)
- The loop assigns higher counts to treatment for DEG_UP genes, lower for DEG_DOWN
- `write.csv` – saves the matrix for later steps

---

# Part 3: Differential Expression with DESeq2

## Page 3.1 – Why DESeq2?

DESeq2 is a widely used R package that:
- Uses a negative binomial model for count data
- Corrects for library size (sequencing depth)
- Shrinks estimates for low-count genes
- Applies multiple-testing correction (padj)

---

## Page 3.2 – Install DESeq2 (First Time Only)

Copy and run:

```r
# Install DESeq2 and other packages (run once)
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "ggplot2", "pheatmap", "RColorBrewer"))
```

---

## Page 3.3 – Run DESeq2 – Full Code with Explanation

```r
# ============================================
# Step: Differential expression analysis
# ============================================

library(DESeq2)

# 1. Load count matrix and metadata
# ---------------------------------
# The count matrix has genes as rows, samples as columns.
# Metadata describes each sample (which condition it belongs to).
# ---------------------------------

counts <- read.csv("04_counts/count_matrix.csv", row.names = 1)
counts <- round(counts)  # DESeq2 needs integer counts

metadata <- data.frame(
  sample = colnames(counts),
  condition = c(rep("Control", 4), rep("Treatment", 4)),
  row.names = colnames(counts)
)

# 2. Create DESeqDataSet
# -----------------------
# design = ~ condition  means: we compare groups defined by "condition"
# ---------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)

# 3. Run DESeq2 analysis
# -----------------------
# This fits the model and tests each gene.
# ---------------------------------

dds <- DESeq(dds)

# 4. Extract results: Treatment vs Control
# ----------------------------------------
# contrast = c("condition", "Treatment", "Control")
#            means: numerator = Treatment, denominator = Control
# log2FC > 0  =>  gene is higher in Treatment
# log2FC < 0  =>  gene is lower in Treatment
# ---------------------------------

res <- results(dds, contrast = c("condition", "Treatment", "Control"))

# 5. Order by adjusted p-value (most significant first)
# -----------------------------------------------------

res_ordered <- res[order(res$padj), ]
res_df <- as.data.frame(res_ordered)
res_df$gene <- rownames(res_df)

# 6. Filter significant DEGs
# ---------------------------
# padj < 0.05       =>  statistically significant
# |log2FC| > 1      =>  at least 2-fold change (up or down)
# ---------------------------------

sig <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)

cat("Number of significant DEGs:", nrow(sig), "\n")
print(head(sig[, c("gene", "log2FoldChange", "padj")]))
```

**Key outputs:**
- `res` – results for all genes
- `log2FoldChange` – fold change (positive = up in Treatment)
- `padj` – adjusted p-value (< 0.05 = significant)
- `sig` – table of significant DEGs

---

# Part 4: Visualizations

## Page 4.1 – Volcano Plot

Shows **fold change** (x-axis) vs **significance** (y-axis). Red points = significant DEGs.

```r
# ============================================
# Volcano plot
# ============================================
# X-axis: log2 Fold Change (left = down, right = up)
# Y-axis: -log10(padj)  (higher = more significant)
# Red points: padj < 0.05 AND |log2FC| > 1
# ============================================

library(ggplot2)

res_df$significant <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
  "Yes", "No"
)
res_df$significant[is.na(res_df$padj)] <- "No"

p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj + 1e-100))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("No" = "gray70", "Yes" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Treatment vs Control",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  )

dir.create("05_deg_results", showWarnings = FALSE)
ggsave("05_deg_results/volcano_plot.pdf", p, width = 6, height = 5)
ggsave("05_deg_results/volcano_plot.png", p, width = 6, height = 5, dpi = 150)
cat("Saved: volcano_plot.pdf and .png\n")
```

---

## Page 4.2 – MA Plot

Shows **mean expression** (x-axis) vs **fold change** (y-axis). Used to check for bias.

```r
# ============================================
# MA plot
# ============================================
# MA = M (log ratio) vs A (mean)
# Helps spot if low-expression genes are spuriously different
# ============================================

p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("No" = "gray70", "Yes" = "red")) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  theme_minimal() +
  labs(
    title = "MA Plot: Treatment vs Control",
    x = "log10(mean expression)",
    y = "log2 Fold Change"
  )

ggsave("05_deg_results/ma_plot.pdf", p_ma, width = 6, height = 5)
cat("Saved: ma_plot.pdf\n")
```

---

## Page 4.3 – Heatmap of Top DEGs

Shows expression of top DEGs across samples. Red = high, blue = low.

```r
# ============================================
# Heatmap of top DEGs
# ============================================
# We use variance-stabilized counts (vst) for better visualization.
# Z-score each gene (row) so colors are comparable.
# ============================================

library(pheatmap)
library(RColorBrewer)

vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation
mat <- assay(vsd)

if (nrow(sig) > 0) {
  top_genes <- rownames(sig)[1:min(20, nrow(sig))]
  mat_heat <- mat[top_genes, ]
  mat_scaled <- t(scale(t(mat_heat)))  # Z-score per gene

  pdf("05_deg_results/heatmap_top_deg.pdf", width = 6, height = 6)
  pheatmap(
    mat_scaled,
    scale = "none",
    color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
    main = "Top DEGs - Z-scored expression",
    border_color = NA,
    fontsize_row = 8
  )
  dev.off()
  cat("Saved: heatmap_top_deg.pdf\n")
}
```

---

# Part 5: Systems Biology

## Page 5.1 – Pathway / Regulation Summary

Groups DEGs into pathways and summarizes up vs down regulation.

```r
# ============================================
# Systems biology: pathway and regulation summary
# ============================================

dir.create("06_systems_biology", showWarnings = FALSE)

if (nrow(sig) > 0) {
  # Regulation: up vs down
  sig$regulation <- ifelse(sig$log2FoldChange > 0, "Upregulated", "Downregulated")
  reg_summary <- as.data.frame(table(sig$regulation))

  p2 <- ggplot(reg_summary, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_col() +
    theme_minimal() +
    labs(title = "DEG Regulation Summary", x = "", y = "Count") +
    scale_fill_brewer(palette = "Set2") +
    theme(legend.position = "none")
  ggsave("06_systems_biology/regulation_summary.pdf", p2, width = 5, height = 4)

  # Save tables
  write.csv(res_df, "05_deg_results/all_genes.csv", row.names = FALSE)
  write.csv(sig, "05_deg_results/significant_deg.csv", row.names = FALSE)
  cat("Saved: all_genes.csv, significant_deg.csv, regulation_summary.pdf\n")
}
```

---

# Part 6: Run Everything in One Go

## Page 6.1 – Complete Script (Copy and Run)

If you prefer to run the entire pipeline at once, use your existing `run_all_in_R.R` or paste this:

```r
# ============================================
# COMPLETE RNA-seq ANALYSIS - Run all steps
# ============================================
# Make sure you are in the tutorial folder (setwd if needed)
# ============================================

# Set your working directory (change path if needed)
setwd("C:/Users/YourName/Desktop/drmahmoodhachim-gif.github.io-main/bioinformatics-tutorial")

# Run the full script
source("run_all_in_R.R")

# Or: if run_all_in_R.R does not exist, the steps above are the full pipeline.
```

---

# Part 7: Quick Reference – Key R Commands

| What you want | R code |
|---------------|--------|
| Install packages | `BiocManager::install(c("DESeq2","ggplot2","pheatmap","RColorBrewer"))` |
| Load DESeq2 | `library(DESeq2)` |
| Read counts | `counts <- read.csv("04_counts/count_matrix.csv", row.names = 1)` |
| Create DDS | `dds <- DESeqDataSetFromMatrix(countData, colData, design = ~ condition)` |
| Run DESeq | `dds <- DESeq(dds)` |
| Get results | `res <- results(dds, contrast = c("condition", "Treatment", "Control"))` |
| Filter DEGs | `sig <- subset(as.data.frame(res), padj < 0.05 & abs(log2FoldChange) > 1)` |

---

# Summary

1. **Design** – Clear groups, enough replicates, balanced.
2. **Count matrix** – Genes × samples, integer counts.
3. **DESeq2** – Statistical testing for differential expression.
4. **Volcano** – Effect size vs significance.
5. **MA plot** – Check for bias.
6. **Heatmap** – Top DEGs across samples.
7. **Pathways** – Biological interpretation.

---

*Good luck with your analysis. For real data, use the same workflow with your own count matrix and metadata.*
