# RNA-seq Pipeline: Student Manual

**DEG Analysis & Systems Biology – No Python Required**

*For 7 PhD students with no programming background, using basic laptops*

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Your Project Assignment](#your-project-assignment)
3. [What You Will Learn](#what-you-will-learn)
4. [Detailed First-Time Setup](#detailed-first-time-setup)
5. [How to Run – Step by Step](#how-to-run--step-by-step)
6. [Understanding Your Results](#understanding-your-results)
7. [Analysis Steps Explained – What Each Step Means](#analysis-steps-explained--what-each-step-means)
8. [Interpreting Results for Your Project](#interpreting-results-for-your-project)
9. [Complete Troubleshooting Guide](#complete-troubleshooting-guide)
10. [Checklist](#checklist)
11. [Folder Structure](#folder-structure)
12. [Next Steps](#next-steps)

---

## Quick Start

**You only need R. No Python. No command line.**

1. **Install R** from https://cran.r-project.org/
2. **Double-click** `RUN_ANALYSIS.bat` (or open `run_all_in_R.R` in RStudio and click **Source**)
3. **View results** in `05_deg_results/volcano_plot.png` and `06_systems_biology/pathway_barplot.png`

**Want to run step-by-step with explanations?** → Open **`RNA-seq_Complete_Manual.md`** for a ChatGPT-style guide: experimental design, each analysis step, and copy-paste R code with explanations.

**Teaching slides:** Open **`slides/index.html`** in your browser for an 11-slide visual overview. Use ← → arrows to navigate.

**BLAST activity:** Open **`BLAST_ACTIVITY.md`** and use **`sequences_for_blast.fa`** to practice sequence similarity search at https://blast.ncbi.nlm.nih.gov (web only, no installation).

---

## Your Project Assignment

| # | Student | Project |
|---|---------|---------|
| 1 | **Mutho** | Epilepsy |
| 2 | **Tuleen** | Cardiac Genetics (Congenital) |
| 3 | **Maryam** | Diabetes Mellitus (DM) |
| 4 | **Zehra** | Breast Cancer – Fibroblasts |
| 5 | **Ayesha** | Breast Cancer – Macrophages |
| 6 | **Ahmed** | Leukemia |
| 7 | **Yara** | Autoimmune Asthma |

**Your folder:** Find your project folder (e.g. `01_MUTHO_Epilepsy`) and read `PROJECT_INFO.txt` for project-specific tips. Save your figures and notes there when done.

---

## What You Will Learn

1. **Count matrix** – Gene expression data from RNA-seq  
2. **Differential expression** – Finding genes that change between conditions  
3. **Visualizations** – Volcano plots, heatmaps, MA plots  
4. **Systems biology** – Pathways and gene regulation  

---

## Detailed First-Time Setup

### Step 1: Download and Extract the Tutorial

1. Download the tutorial folder (ZIP file or folder) from your instructor.
2. **Extract** (if it’s a ZIP): Right-click the ZIP → **Extract All** → choose a simple location like Desktop or Documents.
3. **Important:** Do not put it inside `C:\Program Files`, `OneDrive` sync folders, or very long paths. Use something like:
   - `C:\Users\YourName\Desktop\bioinformatics-rnaseq-tutorial`
   - or `D:\bioinformatics-rnaseq-tutorial`

### Step 2: Install R

1. Go to **https://cran.r-project.org/**
2. Click **"Download R for Windows"** (or **"Download R for macOS"** / **"Download R for Linux"** if you use those).
3. Choose the latest version (e.g. R-4.4.x).
4. Run the installer.
5. **IMPORTANT:** During installation, when you see component options, **check "Add R to the system PATH"** if available. This helps the `.bat` file find R.
6. Use default options for everything else.
7. **Restart your computer** after installing (or at least close and reopen any terminals).

### Step 3: (Optional) Install RStudio

RStudio makes running the script easier and shows errors more clearly.

1. Go to **https://posit.co/download/rstudio-desktop/**
2. Download the free **RStudio Desktop**.
3. Install it (R must be installed first).

### Step 4: Install R Packages (First Run Only)

The first time you run the script, R will likely ask to install these packages:

- DESeq2  
- ggplot2  
- pheatmap  
- RColorBrewer  

**When R asks:** Type **`y`** and press **Enter** to install.

**To install manually** (if something fails), open R or RStudio and paste this, then press Enter:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2","ggplot2","pheatmap","RColorBrewer"))
```

---

## How to Run – Step by Step

### Option 1: Double-Click (Windows)

1. Open File Explorer and go to the `bioinformatics-rnaseq-tutorial` folder.
2. Find **`RUN_ANALYSIS.bat`** (it has a gear or script icon).
3. **Double-click** it.
4. A black window (Command Prompt) will open.
5. You may see:
   - `Found R. Running analysis...` → Good, wait 1–2 minutes.
   - `ERROR: R was not found` → See [Troubleshooting: R Not Found](#1-r-not-found-when-double-clicking-run_analysisbat).
6. When you see **"Analysis complete!"**, the analysis is done.
7. Open the folders:
   - `05_deg_results` → double-click `volcano_plot.png`
   - `06_systems_biology` → double-click `pathway_barplot.png`
8. Press any key in the black window to close it (or close it manually).

### Option 2: RStudio (All Computers – Windows, Mac, Linux)

1. Open **RStudio**.
2. Go to **File → Open File** (or **File → Open**).
3. Browse to the tutorial folder and select **`run_all_in_R.R`**.
4. The script will appear in the editor.
5. Click the **Source** button (top-right of the script panel, or the icon with a document and green play).
6. Watch the Console at the bottom for progress messages.
7. When it says **"DONE! All results saved"**, open `05_deg_results` and `06_systems_biology` to view your figures.

### Option 3: Plain R (No RStudio)

1. Open **R** (the R console, not RStudio).
2. Change to the tutorial folder. In the R console, type (adjust the path if yours is different):

   ```r
   setwd("C:/Users/YourName/Desktop/bioinformatics-rnaseq-tutorial")
   ```

3. Run the script:

   ```r
   source("run_all_in_R.R")
   ```

---

## Understanding Your Results

### What the Script Does (in Order)

1. **Creates synthetic data** – Simulates RNA-seq counts for 50 genes in 8 samples (4 Control, 4 Treatment).
2. **Runs differential expression** – Uses DESeq2 to find genes that change between Control and Treatment.
3. **Creates visualizations** – Volcano plot, heatmap, MA plot.
4. **Systems biology** – Pathway barplot and regulation summary.

### Output Files

| File | Location | What It Is |
|------|----------|------------|
| `volcano_plot.png` | `05_deg_results/` | Volcano plot – genes that change significantly |
| `volcano_plot.pdf` | `05_deg_results/` | Same as PNG, PDF version |
| `ma_plot.pdf` | `05_deg_results/` | MA plot – fold change vs average expression |
| `heatmap_top_deg.pdf` | `05_deg_results/` | Heatmap of top changed genes |
| `significant_deg.csv` | `05_deg_results/` | List of significant genes (open in Excel) |
| `all_genes.csv` | `05_deg_results/` | All genes with statistics |
| `pathway_barplot.png` | `06_systems_biology/` | DEGs by pathway |
| `pathway_barplot.pdf` | `06_systems_biology/` | Same, PDF version |
| `regulation_summary.pdf` | `06_systems_biology/` | Up vs down regulated count |
| `count_matrix.csv` | `04_counts/` | Raw count data (for reference) |

### Volcano Plot

- **X-axis:** log2 Fold Change – how much the gene changed (left = down, right = up).
- **Y-axis:** –log10(adjusted p-value) – statistical significance (higher = more significant).
- **Red points:** Significant genes (padj < 0.05 and |log2FC| > 1).
- **Gray points:** Not significant.

### MA Plot

- **X-axis:** Average expression level.
- **Y-axis:** log2 Fold Change.
- Helps you see if lowly expressed genes are spuriously called as different.

### Heatmap

- **Rows:** Top differentially expressed genes.
- **Columns:** Your 8 samples (CTRL_01–04, TREAT_01–04).
- **Colors:** Expression level (red = high, blue = low after Z-scoring).

---

## Analysis Steps Explained – What Each Step Means

This section explains *what* each analysis step does and *why* it matters for your research.

---

### Step 1: Count Matrix (04_counts/count_matrix.csv)

**What it is**

A table where each row is a gene and each column is a sample. Each cell is the number of RNA-seq reads that mapped to that gene in that sample.

**What it means biologically**

- **Higher counts** = more mRNA from that gene = the gene is more highly expressed in that sample.
- RNA-seq measures how much of each gene’s mRNA is present in the tissue or cells.
- Counts are integer values (1, 2, 3, …) because they count discrete sequencing reads.

**Why it matters**

The count matrix is the main input for differential expression. You compare counts across conditions (e.g. Control vs Treatment) to see which genes change.

**In your project**

- Rows: ~50 genes.
- Columns: 8 samples (4 Control, 4 Treatment).

---

### Step 2: Differential Expression Analysis (DESeq2)

**What it does**

DESeq2 compares gene expression between two conditions (e.g. Control vs Treatment) and tests which genes are significantly different.

**Main outputs**

| Term | Meaning |
|------|---------|
| **log2 Fold Change (log2FC)** | How much expression changes. Positive = higher in Treatment; negative = lower. log2FC = 1 means 2×; log2FC = 2 means 4×. |
| **p-value** | Probability the observed difference is due to chance. Smaller = more likely real. |
| **padj (adjusted p-value)** | p-value corrected for multiple testing (many genes tested). We use padj &lt; 0.05 as significant. |
| **baseMean** | Average expression across all samples. |

**What it means biologically**

- **log2FC &gt; 0** → gene is **upregulated** in Treatment (e.g. disease).
- **log2FC &lt; 0** → gene is **downregulated** in Treatment.
- **padj &lt; 0.05** → we treat the change as statistically significant.

**Why it matters**

DESeq2 accounts for:
- Biological variability between replicates.
- Different sequencing depths across samples.
- Low-count genes that are harder to estimate.

This gives more reliable lists of differentially expressed genes (DEGs).

---

### Step 3: Volcano Plot (volcano_plot.png)

**What it shows**

A scatter plot of all genes with:
- **X-axis:** log2 Fold Change (direction of change).
- **Y-axis:** −log10(padj) (statistical significance).

**How to read it**

- **Upper right:** upregulated and significant (padj &lt; 0.05, log2FC &gt; 1).
- **Upper left:** downregulated and significant (padj &lt; 0.05, log2FC &lt; −1).
- **Dashed horizontal line:** padj = 0.05.
- **Dashed vertical lines:** log2FC = ±1 (2× change).
- **Red points:** genes we call significant DEGs.

**What it means biologically**

- Points far from zero on the x-axis = strong fold changes.
- Points high on the y-axis = strong statistical support.
- Genes in the upper corners are the best candidates for follow-up (e.g. validation, pathway analysis).

**Why we use it**

A volcano plot shows both effect size (fold change) and significance in one view, so you can quickly identify genes that are both large and reliable.

---

### Step 4: MA Plot (ma_plot.pdf)

**What it shows**

- **X-axis:** average expression (log10 scale).
- **Y-axis:** log2 Fold Change.

**How to read it**

- Most points hover around y = 0 (no change).
- Significant DEGs are above or below y = 0.
- We check whether lowly expressed genes (left side) cluster oddly along the y-axis.

**What it means biologically**

- **Good:** DEGs spread across expression levels; no clear bias toward low counts.
- **Warning:** Many DEGs only among lowly expressed genes can indicate batch or technical bias rather than biology.

**Why we use it**

The MA plot helps check if the analysis is biased. DESeq2 shrinks estimates for low-count genes, and the MA plot shows whether that and other assumptions look reasonable.

---

### Step 5: Heatmap (heatmap_top_deg.pdf)

**What it shows**

A grid where:
- **Rows:** top differentially expressed genes.
- **Columns:** samples (CTRL_01–04, TREAT_01–04).
- **Color:** expression level (Z-score: how many standard deviations from the mean).

**How to read it**

- **Red:** higher than average for that gene.
- **Blue:** lower than average for that gene.
- **White:** near average.

**What it means biologically**

- **Vertical patterns:** samples that cluster together have similar expression profiles.
- **Horizontal patterns:** genes that change together may be in the same pathway or regulated by the same factor.
- **Control vs Treatment:** you expect distinct clusters (e.g. CTRL vs TREAT in different blocks).

**Why we use it**

The heatmap summarizes expression of key genes across samples and supports quality checks (e.g. replicates grouping, separation of conditions).

---

### Step 6: Systems Biology – Pathway Barplot & Regulation Summary

**What it shows**

- **Pathway barplot:** how many DEGs fall into each biological pathway (e.g. cell cycle, apoptosis).
- **Regulation summary:** how many DEGs are upregulated vs downregulated.

**What it means biologically**

- **Pathways** are sets of genes that work together (e.g. same metabolic pathway, signaling cascade).
- If many DEGs belong to one pathway, that pathway is likely affected by your condition.
- **Up vs down** tells you whether the condition generally activates or represses gene expression.

**Why it matters**

- Moves you from single genes to biological processes.
- Suggests mechanisms (e.g. “cell cycle dysregulated in disease”).
- Gives context for individual genes you might validate.

**In this tutorial**

We use simple, illustrative pathways. With real data you would use databases (GO, KEGG) and tools like `clusterProfiler` for proper pathway enrichment.

---

### Summary: The Flow of Your Analysis

```
Count matrix (raw data)
    ↓
DESeq2 (statistical testing)
    ↓
Volcano plot, MA plot (visual QC and overview)
    ↓
Heatmap (top DEGs across samples)
    ↓
Pathway analysis (biological meaning)
```

---

## Interpreting Results for Your Project

When you later work with real data in your project, focus on genes relevant to your disease:

| Student | Project | What to Look For |
|---------|---------|------------------|
| **Mutho** | Epilepsy | Ion channels, synaptic genes, seizure pathways |
| **Tuleen** | Cardiac Genetics | Heart development genes (NKX2-5, TBX5, GATA4) |
| **Maryam** | DM | Insulin signaling, beta-cells, diabetes inflammation |
| **Zehra** | Breast Cancer Fibroblasts | CAF markers, ECM, tumor stroma |
| **Ayesha** | Breast Cancer Macrophages | TAM markers, M1/M2, immune suppression |
| **Ahmed** | Leukemia | Hematopoiesis, leukemia drivers, cell cycle |
| **Yara** | Autoimmune Asthma | Th2/Th17, airway remodeling, allergic inflammation |

---

## Complete Troubleshooting Guide

### 1. R Not Found When Double-Clicking RUN_ANALYSIS.bat

**Symptoms:** Window says "ERROR: R was not found" or "R is not installed".

**Solutions (try in order):**

- **A. Reinstall R and add to PATH**
  1. Uninstall R (Control Panel → Uninstall).
  2. Download R again from https://cran.r-project.org/.
  3. During install, look for **"Add R to the system PATH"** or **"Modify the system PATH"** and **check it**.
  4. Restart your computer.
  5. Try the .bat file again.

- **B. Use RStudio instead**
  1. Install R from https://cran.r-project.org/.
  2. Install RStudio from https://posit.co/download/rstudio-desktop/.
  3. Open RStudio → File → Open → select `run_all_in_R.R`.
  4. Click **Source**. (You can ignore the .bat file.)

- **C. Move the tutorial folder**
  - Some folder locations (e.g. OneDrive, network drives) can cause issues.
  - Copy the tutorial to `C:\Users\YourName\Desktop\bioinformatics-rnaseq-tutorial` and try again.

---

### 2. R Asks to Install Packages – "Do you want to install from sources?"

**Symptoms:** R says "Do you want to install from sources?" or "These packages need to be installed".

**Solution:** Type **`n`** (no) and press Enter. R will usually install from a binary. If it still fails, type **`y`** and wait (it may take 5–10 minutes).

---

### 3. "there is no package called 'DESeq2'" (or ggplot2, pheatmap, RColorBrewer)

**Symptoms:** Error like `there is no package called 'DESeq2'` or `could not find package "DESeq2"`.

**Solution:** Open R or RStudio and run:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2","ggplot2","pheatmap","RColorBrewer"))
```

Wait for installation to finish (may take a few minutes). Then run the script again.

---

### 4. BiocManager or Package Installation Fails

**Symptoms:** Error when running `BiocManager::install(...)`, e.g. "installation of package failed" or "had non-zero exit status".

**Solutions:**

- **A. Try updating R**
  - Install the latest R from https://cran.r-project.org/.
  - Older R versions sometimes fail with newer packages.

- **B. Install one package at a time**
  ```r
  install.packages("ggplot2")
  install.packages("pheatmap")
  install.packages("RColorBrewer")
  BiocManager::install("DESeq2")
  ```

- **C. Install dependencies first**
  ```r
  install.packages(c("Matrix", "S4Vectors", "IRanges", "GenomicRanges"))
  BiocManager::install("DESeq2")
  ```

- **D. If you have antivirus/firewall**
  - Temporarily disable it or add an exception for R, then try installing again.

---

### 5. Window Closes Too Fast / I Can't See the Output

**Symptoms:** Black window appears and disappears quickly; you can't read any messages.

**Solution:** Run via RStudio instead:

1. Open RStudio.
2. File → Open → `run_all_in_R.R`.
3. Click **Source**.
4. All messages appear in the Console and stay visible.

---

### 6. "Access is denied" or "Permission denied"

**Symptoms:** Error about access or permission when running the script.

**Solutions:**

- **A. Avoid system folders**
  - Don’t put the tutorial in `C:\Program Files`, `C:\Windows`, or `C:\system32`.

- **B. Use a simple path**
  - Move the folder to `C:\Users\YourName\Desktop\` or `C:\Users\YourName\Documents\`.

- **C. Run as normal user**
  - Right-click the folder → Properties → ensure it’s not read-only.
  - If on a shared PC, copy the folder to your personal Documents/Desktop.

---

### 7. Can't Find My Results / No 05_deg_results Folder

**Symptoms:** The script seemed to run, but you don’t see `05_deg_results` or `06_systems_biology`.

**Solutions:**

- **A. Check folder location**
  - Results are created **inside** the `bioinformatics-rnaseq-tutorial` folder.
  - Make sure you’re looking in the correct folder (the one that contains `run_all_in_R.R`).

- **B. Check if the script actually finished**
  - In RStudio: look at the Console for "DONE! All results saved".
  - If you see an error before that, the script stopped early.

- **C. Run again**
  - Close R/RStudio, reopen, and run the script again. The folders will be created.

---

### 8. Script Stops with an Error (Red Text in RStudio)

**Symptoms:** Red error message in the RStudio Console; script stops.

**What to do:**

1. **Read the error message** – it often says which package or function failed.
2. **Common errors and fixes:**
   - `could not find function "ggplot"` → Install ggplot2: `install.packages("ggplot2")`
   - `object 'xxx' not found` → Usually means a step failed earlier. Run the whole script from the top (Source).
   - `cannot open file` → The working directory is wrong. In R: `setwd("C:/path/to/bioinformatics-rnaseq-tutorial")` then run again.
   - `subscript out of bounds` or `replacement has X rows, data has Y` → Unlikely with this script; try running it again.

3. **Copy the full error** and send it to your instructor if you’re still stuck.

---

### 9. Analysis Runs But Figures Look Wrong / Empty

**Symptoms:** Folders exist but PNGs are blank, or plots look strange.

**Solutions:**

- **A. Antivirus blocking**
  - Some antivirus software blocks R from writing files.
  - Add an exception for R or the tutorial folder, or temporarily disable it.

- **B. Run again**
  - Sometimes the first run fails partway through.
  - Close R/RStudio and run the script again.

---

### 10. "R version 4.x is required" or Package Version Errors

**Symptoms:** Error like "this package was built for R version 4.x" or "R version too old".

**Solution:** Install the latest R from https://cran.r-project.org/ (e.g. R-4.4.x or newer).

---

### 11. Slow / Freezing / Taking Very Long

**Symptoms:** Script runs for a long time or seems to hang.

**Note:** On a normal laptop, this script usually finishes in 1–2 minutes.

**Solutions:**

- Close other programs (browsers, Word, etc.) to free memory.
- If you have very little RAM (e.g. 4 GB), close everything except R/RStudio.
- Wait up to 5 minutes; package loading can be slow the first time.

---

### 12. Mac or Linux – .bat File Doesn't Work

**Symptoms:** On Mac or Linux, double-clicking `RUN_ANALYSIS.bat` does nothing.

**Solution:** Use the shell script or RStudio:

```bash
cd /path/to/bioinformatics-rnaseq-tutorial
bash RUN_ANALYSIS.sh
```

Or open `run_all_in_R.R` in RStudio and click **Source**.

---

### 13. "The system cannot find the path specified"

**Symptoms:** Error about path not found when running the .bat file.

**Solutions:**

- Don’t use paths with special characters (é, ñ, 中文, etc.).
- Don’t use very long paths – keep the folder name short.
- Move the folder to `C:\Users\YourName\Desktop\bioinf` and try again.

---

### Quick Troubleshooting Reference

| Problem | Quick Fix |
|---------|-----------|
| R not found | Use RStudio: open `run_all_in_R.R` → Source |
| Package missing | In R: `BiocManager::install(c("DESeq2","ggplot2","pheatmap","RColorBrewer"))` |
| Window closes fast | Use RStudio instead of .bat |
| No results folder | Check you're in the right folder; run script again |
| Access denied | Move folder to Desktop or Documents |
| Script error | Read error message; install missing package; run again |

---

## Checklist

- [ ] Downloaded and extracted the tutorial folder  
- [ ] Installed R (and optionally RStudio)  
- [ ] Opened `PROJECT_INFO.txt` in my project folder  
- [ ] Ran the analysis (double-click `RUN_ANALYSIS.bat` OR Source in RStudio)  
- [ ] Opened `volcano_plot.png` and `pathway_barplot.png`  
- [ ] Copied my figures to my project folder  

---

## Folder Structure

```
bioinformatics-rnaseq-tutorial/
├── START_HERE.txt
├── RUN_ANALYSIS.bat         ← Double-click (Windows)
├── RUN_ANALYSIS.sh          ← Use on Mac/Linux
├── run_all_in_R.R           ← Or open in RStudio and Source
├── README.md
├── STUDENT_MANUAL.md        ← This file
├── 04_counts/               (created when you run)
├── 05_deg_results/          (volcano, heatmap, DEG list)
├── 06_systems_biology/      (pathway plots)
├── 01_MUTHO_Epilepsy/
├── 02_TULEEN_CardiacGeneticsCongenital/
├── 03_MARYAM_DM/
├── 04_ZEHRA_BreastCancerFibroblasts/
├── 05_AYESHA_BreastCancerMacrophages/
├── 06_AHMED_Leukemia/
└── 07_YARA_AutoimmuneAsthma/
```

---

## Next Steps

1. Use this workflow on real data from your project.  
2. Learn GO/KEGG enrichment with `clusterProfiler` in R for real gene data.  

---

*Good luck with your analysis. If you get stuck, check the troubleshooting section or ask your instructor.*
