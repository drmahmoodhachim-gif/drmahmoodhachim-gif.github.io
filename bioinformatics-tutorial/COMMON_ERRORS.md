# Common Errors & Fixes

If you get an error, find it in this list and follow the fix.

---

## 1. "R was not found" or "R is not installed"

**What you see:** Black window says "ERROR: R was not found"

**Fix:**
1. Install R from https://cran.r-project.org/
2. **Or** use RStudio: Open RStudio → File → Open → `run_all_in_R.R` → Click **Source**

---

## 2. "there is no package called 'DESeq2'" (or ggplot2, pheatmap, RColorBrewer)

**What you see:** Red error in R: `there is no package called 'DESeq2'`

**Fix:** In R or RStudio, run this first:
```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2","ggplot2","pheatmap","RColorBrewer"))
```
Then run the analysis again.

---

## 3. "cannot change working directory" or "cannot open file '04_counts/count_matrix.csv'"

**What you see:** Error about file not found or wrong directory

**Fix:** You must run the script from the **bioinformatics-tutorial** folder (the one containing run_all_in_R.R), **not** from 01_MUTHO_Epilepsy or other student folders.

- **If using .bat:** Double-click `RUN_ANALYSIS.bat` from inside the bioinformatics-tutorial folder.
- **If using RStudio:** Open `run_all_in_R.R` from bioinformatics-tutorial. Before Source, run:
  ```r
  setwd("C:/Users/YOUR_USERNAME/Desktop/drmahmoodhachim-gif.github.io-main/bioinformatics-tutorial")
  ```
  (Change the path to where YOUR bioinformatics-tutorial folder is. Use forward slashes.)

---

## 4. "Error in library(DESeq2) : package 'DESeq2' is not available"

**What you see:** DESeq2 won't load even after trying to install

**Fix:** Your R version may be too old. Install the latest R from https://cran.r-project.org/ (e.g. R-4.4.x). Then:
```r
BiocManager::install("DESeq2")
```

---

## 5. R asks "Do you want to install from sources?" or "Update all/some/none?"

**What you see:** Prompt during package installation

**Fix:** Type **n** (no) and press Enter. R will try to install a pre-built version. If it still fails, type **y** and wait (may take 5–10 minutes).

---

## 6. "attempt to set 'rownames' on an object with no dimensions"

**What you see:** Error at Step 1, when creating synthetic data

**Fix:** This was a bug in older versions of the tutorial – now fixed. Re-download the tutorial from the website to get the latest version. If you already have the latest, the script should work.

---

## 7. "less than 'nsub' rows" (vst error)

**What you see:** Error in vst(): "less than 'nsub' rows, it is recommended to use varianceStabilizingTransformation"

**Fix:** This was a bug with small datasets (50 genes) – now fixed in the latest version. Re-download the tutorial. The script now uses a fallback for small gene sets.

---

## 8. "subscript out of bounds" or "replacement has X rows, data has Y"

**What you see:** Red error with those words

**Fix:** Usually means the script ran partially. Close R/RStudio completely, reopen, and run the script again from the beginning (Source).

---

## 9. Window closes immediately / I can't read the error

**What you see:** Black window flashes and closes

**Fix:** Use RStudio instead:
1. Open RStudio
2. File → Open → `run_all_in_R.R`
3. Click **Source**
4. All messages (and errors) will stay visible in the Console

---

## 10. "Access is denied" or "Permission denied"

**What you see:** Error when writing files

**Fix:**
- Don't run from `C:\Program Files` or `C:\Windows`
- Copy the tutorial folder to your **Desktop** or **Documents**
- Make sure the folder is not read-only (right-click → Properties → uncheck Read-only)

---

## 11. Antivirus blocking or files not created

**What you see:** Script says "DONE" but no PDF/PNG files, or antivirus message

**Fix:** Add an exception for R (or the tutorial folder) in your antivirus. Or temporarily disable antivirus while running, then turn it back on.

---

## 12. "object 'sig' not found" or "object 'res_df' not found"

**What you see:** Error about object not found

**Fix:** You ran part of the script but not all. Run the **entire** script from the top (click Source, or run `source("run_all_in_R.R")`).

---

---

## Mac-Specific Issues

### "RUN_ANALYSIS.bat" does nothing on Mac

**Cause:** `.bat` files are for Windows only.

**Fix:** Use **RStudio** (File → Open → run_all_in_R.R → Source) or see **MAC_INSTRUCTIONS.md** for the Terminal method.

---

### Mac: "Rscript: command not found" in Terminal

**Fix:** Use **RStudio** instead – it finds R automatically. Or add R to your PATH (R installs to `/Library/Frameworks/R.framework/`).

---

### Mac: "Permission denied" when running RUN_ANALYSIS.sh

**Fix:** In Terminal:
```bash
chmod +x RUN_ANALYSIS.sh
./RUN_ANALYSIS.sh
```

---

### Mac: Double-clicking .sh file opens it in TextEdit

**Fix:** Use **RStudio** or run from Terminal: `bash RUN_ANALYSIS.sh` (see MAC_INSTRUCTIONS.md).

---

### Mac: "Setting LC_CTYPE failed" / "non-UTF8 locale" warnings

**What you see:** Warnings when R starts: "Setting LC_CTYPE failed...", "You're using a non-UTF8 locale..."

**Fix:** The script usually still runs. If you get errors, run this in R first:
```r
Sys.setenv(LANG = "en_US.UTF-8")
Sys.setenv(LC_ALL = "en_US.UTF-8")
```
Then run the script. Or see MAC_INSTRUCTIONS.md for the Terminal method.

---

## 13. "What is the second part? / DEGs in the individual projects?"

**What you're wondering:** The instructions mention "interpret your DEGs in the context of your project" – does that mean a separate analysis for my project folder?

**Answer:** No. There is **only one pipeline** – you run it once and get volcano_plot.png, pathway_barplot.png, and significant_deg.csv. The "second part" is **interpretation**, not a second run.

**What to do:**
1. Run the pipeline (Part 1) → you get the 3 outputs.
2. **Interpret** those results in the context of your project theme (Part 2):
   - Open your **PROJECT_INFO.txt** (e.g. 01_MUTHO_Epilepsy) – it lists what to look for (e.g. ion channels for epilepsy).
   - Look at your volcano plot, pathway barplot, and significant_deg.csv.
   - Write 2–3 sentences: "I found X DEGs. Pathway Y is enriched, which in my project (epilepsy/diabetes/etc.) could relate to …"

**Example (Epilepsy):** "I found 10 DEGs. The pathway barplot shows Cell Cycle is enriched – in epilepsy this could relate to neuronal changes after seizures."

The genes (GENE_001, GENE_002, etc.) are synthetic – the *interpretation* in your project context is what makes it "individual" to your assignment.

---

## 14. Submission form won't let me upload

**What you see:** The Submit button goes nowhere, or the form doesn't accept files.

**Fix:** Email your instructor directly. Attach volcano_plot.png, pathway_barplot.png, and significant_deg.csv. In the email body, include your name, project, number of DEGs, and 2–3 sentences about your findings. This is acceptable.

---

## Still stuck?

1. **Copy the full error message** (select the red text, Cmd+C on Mac)
2. **Take a screenshot** of the R/RStudio window
3. **Send both** to your instructor with a note about what you did (e.g. "I double-clicked RUN_ANALYSIS.bat")
