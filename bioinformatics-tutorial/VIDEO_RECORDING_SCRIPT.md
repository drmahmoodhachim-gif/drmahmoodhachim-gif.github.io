# Video Recording Script for RNA-seq Tutorial

Use this script to record a screen-capture walkthrough for your students. Suggested tools: **OBS Studio** (free), **Loom** (free tier), or **Windows Game Bar** (Win+G).

---

## Video 1: Overview + How to Download and Run (10–12 min)

### Intro (0:00–1:00)
> "Hi everyone. This short video will walk you through the RNA-seq tutorial for your PhD projects. You'll learn how to run a differential gene expression analysis from start to finish, with no Python and no command line. Let's get started."

### Show the Website (1:00–2:30)
> "First, open the tutorial website. Go to drmahmoodhachim-gif.github.io."
- Open https://drmahmoodhachim-gif.github.io in a browser
> "Here you see the main page with Quick Start steps, project assignments for each of you, and documentation links."

> "To download, click the blue 'Download ZIP' button. No login is needed."
- Click Download ZIP

### Extract and Locate (2:30–4:00)
> "Once the ZIP is downloaded, right-click it and choose 'Extract All'. Pick a simple location like your Desktop."
- Extract the ZIP
> "After extracting, you'll see a folder. Open it, then open the **bioinformatics-tutorial** folder inside. This is where you'll work."

### Install R (4:00–6:00)
> "You need R installed. Go to cran.r-project.org. Click 'Download R for Windows' — or for Mac or Linux if that's what you use."
- Go to cran.r-project.org
> "Download and run the installer. During installation, if you see an option to 'Add R to the system PATH', check it. Then finish the installation."

### Run the Analysis (6:00–9:00)
> "Back in the bioinformatics-tutorial folder, find the file RUN_ANALYSIS.bat. Double-click it."
- Double-click RUN_ANALYSIS.bat
> "A black window will open. R will run the analysis. The first time, it may ask to install some packages. Type 'y' and press Enter if it does."
> "Wait about one to two minutes. When you see 'Analysis complete!', you're done."

### Show Results (9:00–11:00)
> "Open the folder 05_deg_results. Double-click volcano_plot.png. This is your volcano plot — red points are genes that change significantly between conditions."
- Open volcano_plot.png
> "Also open 06_systems_biology and look at pathway_barplot.png. This shows which biological pathways your significant genes belong to."

### Close (11:00–12:00)
> "For more details, read START_HERE.txt and RNA-seq_Complete_Manual.md in the tutorial folder. Good luck with your analysis."

---

## Video 2: Understanding the Results (8–10 min)

### Intro (0:00–0:30)
> "In this video we'll explain what each result means and how to interpret your figures."

### Volcano Plot (0:30–2:30)
> "The volcano plot has two axes. On the x-axis is log2 Fold Change — how much each gene changed. Positive values mean the gene is higher in treatment; negative means lower."
> "On the y-axis is minus log10 of the adjusted p-value — that's statistical significance. Higher means more significant."
> "Red points are genes we call significant: they have padj less than 0.05 and an absolute fold change greater than 1. Those are your differentially expressed genes, or DEGs."

### MA Plot (2:30–4:00)
> "The MA plot shows average expression on the x-axis and fold change on the y-axis. Most genes cluster around zero on the y-axis — no change. The MA plot helps us check that we're not over-calling changes in lowly expressed genes."

### Heatmap (4:00–5:30)
> "The heatmap shows our top DEGs across samples. Red means higher expression, blue means lower. Rows are genes, columns are samples. We expect Control samples to cluster together and Treatment samples to cluster together."

### Systems Biology (5:30–7:00)
> "The pathway barplot groups DEGs by biological pathway — for example, cell cycle or metabolism. The regulation summary shows how many genes are upregulated versus downregulated."

### Wrap-up (7:00–8:00)
> "For step-by-step R code and more detail, see RNA-seq_Complete_Manual.md. You can copy and paste the code and run each part yourself."

---

## Suggested Video Links to Add (Optional)

You can also share these existing resources with students:

| Topic | Suggested Video | Link (search) |
|-------|-----------------|---------------|
| RNA-seq overview | StatQuest – RNA-seq | YouTube: "StatQuest RNA seq" |
| DESeq2 | StatQuest – DESeq2 | YouTube: "StatQuest DESeq2" |
| Volcano plots | General bioinformatics | YouTube: "volcano plot RNA seq" |

---

## Recording Tips

1. **Resolution**: 1920×1080 if possible.
2. **Speed**: Speak clearly; students can pause or rewatch.
3. **Cursor**: Ensure the mouse cursor is visible.
4. **Audio**: Use a quiet room and a decent microphone.
5. **Length**: 10–15 minutes per video is usually enough.
