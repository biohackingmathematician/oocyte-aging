#!/bin/bash
# Run All Visualizations Script
# This script generates all visualization files

echo "Running all visualizations..."
echo ""

# Change to script directory
cd "$(dirname "$0")"

echo "Generating forum EDA visualizations (PCA by stage & health score)..."
python3 create_forum_eda_visualizations.py
echo ""

echo "Generating raw datasets EDA..."
python3 create_raw_datasets_eda.py
echo ""

echo "Generating forum raw EDA visualizations..."
python3 create_forum_raw_eda.py
echo ""

echo "Generating complete results summary..."
python3 create_complete_results_summary.py
echo ""

echo "All visualizations generated."
echo ""
echo "Generated files are in ../visualizations/"
ls -1 ../visualizations/*.png 2>/dev/null | wc -l | xargs echo "Total visualizations:"

