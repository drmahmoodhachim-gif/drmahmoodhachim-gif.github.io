# RNA-seq Pipeline Tutorial for PhD Students

**DEG analysis & systems biology – no Python required**

Designed for **7 PhD students** with **no bioinformatics background** and **basic laptops**.  
**Students only need R** – no Python, no command line experience needed.

## Simplest way to run (recommended for students)

1. Install R: https://cran.r-project.org/
2. Double-click **`RUN_ANALYSIS.bat`** (Windows) or run **`bash RUN_ANALYSIS.sh`** (Mac/Linux)
3. Open `05_deg_results/volcano_plot.png` to see results

Or: Open `run_all_in_R.R` in RStudio and click **Source**.

## Student Project Folders & Manual

| Student | Project | Folder |
|---------|---------|--------|
| Mutho | Epilepsy | `01_MUTHO_Epilepsy/` |
| Tuleen | Cardiac Genetics (Congenital) | `02_TULEEN_CardiacGeneticsCongenital/` |
| Maryam | Diabetes Mellitus | `03_MARYAM_DM/` |
| Zehra | Breast Cancer (Fibroblasts) | `04_ZEHRA_BreastCancerFibroblasts/` |
| Ayesha | Breast Cancer (Macrophages) | `05_AYESHA_BreastCancerMacrophages/` |
| Ahmed | Leukemia | `06_AHMED_Leukemia/` |
| Yara | Autoimmune Asthma | `07_YARA_AutoimmuneAsthma/` |

**→ Each student: Read `START_HERE.txt` first, then `STUDENT_MANUAL.md`.**

---

## Quick Start (R only – no Python)

### 1. Install R (one-time)
Download from https://cran.r-project.org/ and run the installer.

### 2. Run the analysis
**Windows:** Double-click `RUN_ANALYSIS.bat`  
**Mac/Linux:** `bash RUN_ANALYSIS.sh`  
**Or:** Open `run_all_in_R.R` in RStudio and click **Source**

### 3. View results
Open `05_deg_results/volcano_plot.png` and `06_systems_biology/pathway_barplot.png`.

---

## What the script does

1. Creates synthetic RNA-seq count data (8 samples: 4 Control, 4 Treatment)
2. Runs differential expression analysis (DESeq2)
3. Generates volcano plot, heatmap, MA plot
4. Creates systems biology visualizations (pathway barplot, regulation summary)

**Output folders** (created when you run): `04_counts/`, `05_deg_results/`, `06_systems_biology/`

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| R not found | Install R from cran.r-project.org, or use RStudio and open run_all_in_R.R → Source |
| R package errors | In R: `BiocManager::install(c("DESeq2","ggplot2","pheatmap","RColorBrewer"))` |
| Out of memory | Close other apps; data is small |

---

## License

Free for educational use.
