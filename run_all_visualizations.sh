#!/bin/bash
# Run All Visualizations Script
# This script generates all visualization files

echo "=========================================="
echo "RUNNING ALL VISUALIZATIONS"
echo "=========================================="
echo ""

# Change to script directory
cd "$(dirname "$0")"

echo "[1/3] Generating forum EDA visualizations (PCA by stage & health score)..."
python3 create_forum_eda_visualizations.py
echo ""

echo "[2/3] Generating improved forum visualizations..."
python3 review_and_fix_visualizations.py
echo ""

echo "[3/3] Fixing existing visualizations (risk, gplvm, complete summary)..."
python3 fix_existing_visualizations.py
echo ""

echo "=========================================="
echo "✓ ALL VISUALIZATIONS GENERATED!"
echo "=========================================="
echo ""
echo "Generated files:"
ls -1 *fixed*.png forum*.png 2>/dev/null | sort | while read f; do
    size=$(ls -lh "$f" 2>/dev/null | awk '{print $5}')
    echo "  ✓ $f ($size)"
done
echo ""
echo "All visualizations are ready for your project forum!"

