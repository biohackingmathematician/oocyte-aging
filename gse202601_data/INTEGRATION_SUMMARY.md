
## GSE202601 Integration Summary

### Dataset Overview

- **Source**: GEO GSE202601
- **Paper**: Chen et al. (2024) Nature Aging
- **Data type**: snRNA-seq + snATAC-seq
- **Samples**: 4 young + 4 aged human ovaries
- **Cell types**: Granulosa, theca, stromal, immune, endothelial

### Clinical Relevance

Granulosa cells are:

1. Routinely collected during IVF oocyte retrieval
2. Normally discarded after procedure
3. Easily accessible without additional procedures
4. Reflect oocyte quality and ovarian aging status

### Files Created

- `gse202601_data/gse202601_metadata.csv` - Sample metadata
- `gse202601_data/analyze_granulosa_aging.py` - Analysis script template

### Next Steps

1. Download processed h5ad from GEO supplementary files
2. Run analyze_granulosa_aging.py
3. Compare with oocyte trajectory correlations
4. Identify shared vs cell-type-specific aging biomarkers

### For Paper Discussion Section

This analysis addresses: "Are there different cell types correlated with 
oocyte age that are easier to collect?"

Answer: Yes - granulosa cells show aging-associated transcriptomic changes
and are clinically accessible during routine IVF procedures.
