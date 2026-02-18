#!/bin/bash
# RNA-seq Tutorial - Run analysis (Mac/Linux)
# No Python needed - R only!
# Usage: Double-click this file, or run: bash RUN_ANALYSIS.sh

cd "$(dirname "$0")"
echo ""
echo "========================================"
echo "  RNA-seq Pipeline for Students"
echo "  (No Python needed - R only!)"
echo "========================================"
echo ""

if ! command -v Rscript &> /dev/null; then
    echo "ERROR: R is not installed."
    echo "Install R from: https://cran.r-project.org/"
    exit 1
fi

echo "Running analysis... (1-2 minutes)"
echo ""
Rscript run_all_in_R.R
echo ""
echo "Done! Check 05_deg_results/ and 06_systems_biology/"
