# Mac Instructions – RNA-seq Tutorial

**If you have a Mac, follow these steps.**

⚠️ **Do not use RUN_ANALYSIS.bat** – it is for Windows only and will not work on Mac.

---

## Quick Start (Mac)

1. **Install R** from https://cran.r-project.org/ → click **"Download R for macOS"**
2. **Install RStudio** from https://posit.co/download/rstudio-desktop/ (recommended)
3. **Open RStudio** → File → Open → select `run_all_in_R.R` from the bioinformatics-tutorial folder
4. **Click Source** (or Cmd+Shift+Enter)
5. Wait 1–2 minutes, then check `05_deg_results/` and `06_systems_biology/`

---

## Option 1: RStudio (Easiest for Mac)

1. Open **RStudio**
2. **File** → **Open File** → go to your `bioinformatics-tutorial` folder → select **run_all_in_R.R**
3. Click the **Source** button (top-right of the script panel)
4. Watch the Console for progress. When you see "DONE!", open the results folders.

**If you get "cannot open file" or "wrong directory":**  
In the R Console, run (change the path to match your folder):

```r
setwd("/Users/YOUR_USERNAME/Desktop/drmahmoodhachim-gif.github.io-main/bioinformatics-tutorial")
```

Then click **Source** again.

---

## Option 2: Terminal + Shell Script

1. Open **Terminal** (Applications → Utilities → Terminal)
2. Go to the tutorial folder:
   ```bash
   cd ~/Desktop/drmahmoodhachim-gif.github.io-main/bioinformatics-tutorial
   ```
   (Adjust the path if your folder is elsewhere – e.g. in Downloads)
3. Make the script executable (first time only):
   ```bash
   chmod +x RUN_ANALYSIS.sh
   ```
4. Run it:
   ```bash
   ./RUN_ANALYSIS.sh
   ```
   Or: `bash RUN_ANALYSIS.sh`

---

## Common Mac Issues

### Locale warnings: "Setting LC_CTYPE failed" / "non-UTF8 locale"

**What you see:** Warnings in R like "Setting LC_CTYPE failed, using 'C'" and "You're using a non-UTF8 locale, therefore only ASCII characters will work."

**Fix – run this in R before clicking Source:**
```r
Sys.setenv(LANG = "en_US.UTF-8")
Sys.setenv(LC_ALL = "en_US.UTF-8")
```

Or in **Terminal** (before opening R):
```bash
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
```

The script will usually still run despite these warnings. If you get errors with special characters, the above fixes help.

---

### "Rscript: command not found"

**Cause:** R is installed but not in the Terminal PATH.

**Fix:** Use **RStudio** instead (Option 1 above). RStudio finds R automatically. You don't need the Terminal.

---

### "Permission denied" when running RUN_ANALYSIS.sh

**Fix:** Run:
```bash
chmod +x RUN_ANALYSIS.sh
./RUN_ANALYSIS.sh
```

---

### Double-clicking RUN_ANALYSIS.sh opens it in TextEdit

**Cause:** Mac opens .sh files in a text editor by default.

**Fix:** Use Terminal (Option 2) or RStudio (Option 1). Don't double-click the .sh file.

---

### "run_all_in_R.R" cannot be opened because it is from an unidentified developer

**Cause:** Mac security (Gatekeeper) blocking files from the internet.

**Fix:**  
- Right-click the file → **Open** → **Open** again to confirm.  
- Or: System Settings → Privacy & Security → allow the app.

---

### "cannot open file '04_counts/count_matrix.csv'"

**Cause:** Working directory is wrong (e.g. you opened R from the wrong folder).

**Fix:** In RStudio, before clicking Source, run:
```r
setwd("/Users/YOUR_USERNAME/Desktop/drmahmoodhachim-gif.github.io-main/bioinformatics-tutorial")
```
Replace `YOUR_USERNAME` with your Mac username. Use the path where your folder actually is (e.g. Downloads).

---

## Path on Mac vs Windows

| Windows           | Mac                          |
|-------------------|------------------------------|
| `C:\Users\Name\Desktop\` | `/Users/Name/Desktop/` |
| Backslash `\`     | Forward slash `/`            |

In R, always use **forward slashes** `/` in paths, even on Windows.

---

## Summary

**Recommended for Mac:** Use **RStudio** + **Source** on `run_all_in_R.R`. No Terminal or shell script needed.
