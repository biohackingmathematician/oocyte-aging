#!/bin/bash
# Run All Visualizations Script
# This script generates all visualization files

echo "=========================================="
echo "RUNNING ALL VISUALIZATIONS"
echo "=========================================="
echo ""

# Change to script directory
cd "$(dirname "$0")"

echo "[1/4] Generating forum EDA visualizations (PCA by stage & health score)..."
python3 create_forum_eda_visualizations.py
echo ""

echo "[2/4] Generating raw datasets EDA..."
python3 create_raw_datasets_eda.py
echo ""

echo "[3/4] Generating forum raw EDA visualizations..."
python3 create_forum_raw_eda.py
echo ""

echo "[4/4] Generating complete results summary..."
python3 create_complete_results_summary.py
echo ""

echo "=========================================="
echo "✓ ALL VISUALIZATIONS GENERATED!"
echo "=========================================="
echo ""
echo "Generated files are in ../visualizations/"
ls -1 ../visualizations/*.png 2>/dev/null | wc -l | xargs echo "Total visualizations:"
echo ""
echo "All visualizations are ready!"

