#!/bin/bash

# Script to run all Snakemake files in the 06-homolog-filter pipeline
# The files are run in numerical order based on their prefix

set -e  # Exit on error

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "========================================="
echo "Running all Snakemake workflows in order"
echo "========================================="
echo ""

# Step 1: Filter KEGG IDs (k01)
echo "[1/6] Running 00-filter-kegg-k01.smk..."
pixi run snakemake -s 00-filter-kegg-k01.smk -c90 -F
echo "✓ Completed 00-filter-kegg-k01.smk"
echo ""

# Step 2: Find domain and filter IDs
echo "[2/6] Running 01-find-domain-filter-ids.smk..."
pixi run snakemake -s 01-find-domain-filter-ids.smk -c90 -F
echo "✓ Completed 01-find-domain-filter-ids.smk"
echo ""

# Step 3: Filter by domain presence (f00/f01)
echo "[3/6] Running 02-filter-domain-f00-f01.smk..."
pixi run snakemake -s 02-filter-domain-f00-f01.smk -c90 -F
echo "✓ Completed 02-filter-domain-f00-f01.smk"
echo ""

# Step 4: Filter eukaryotic copies (e01)
echo "[4/6] Running 03-filter-euk-copies-e01.smk..."
pixi run snakemake -s 03-filter-euk-copies-e01.smk -c90 -F
echo "✓ Completed 03-filter-euk-copies-e01.smk"
echo ""

# Step 5: CD-HIT clustering (c01)
echo "[5/6] Running 04-cdhit-c01.smk..."
pixi run snakemake -s 04-cdhit-c01.smk -c90 -F
echo "✓ Completed 04-cdhit-c01.smk"
echo ""

# Step 6: Manual filter m01 
echo "[6/6] Running 05-filter-manual-m01.smk..."
pixi run snakemake -s 05-filter-manual-m01.smk -c90 -F
echo "✓ Completed 05-filter-manual-m01.smk"
echo ""

# Step 7: Manual filter m02 
echo "[7/7] Running 06-filter-manual-m02.smk..."
pixi run snakemake -s 06-filter-manual-m02.smk -c90 -F
echo "✓ Completed 06-filter-manual-m02.smk"
echo ""

echo "========================================="
echo "All Snakemake workflows completed!"
echo "========================================="
