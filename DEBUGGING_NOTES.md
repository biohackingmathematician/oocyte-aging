# Debugging Notes for Upgrade Sections

## Issues Found and Fixed

### 1. Package Installation
- **Issue**: Packages may not be installed
- **Fix**: Added comprehensive package installation cell (Cell 65)
- **Action**: Run Cell 65 first to install missing packages

### 2. Missing adata Object
- **Issue**: Sections assume `adata` exists
- **Fix**: Added pre-flight checks cell (Cell 66) that validates `adata` exists
- **Action**: Ensure data loading cells have been run before upgrade sections

### 3. Import Errors
- **Issue**: scvi, gpflow, tensorflow may not be available
- **Fix**: All sections have try-except blocks for graceful degradation
- **Action**: Install missing packages: `pip install scvi-tools gpflow tensorflow`

## Running Order

1. **Cell 65**: Install required packages
2. **Cell 66**: Pre-flight checks (validates adata and packages)
3. **Cell 67-68**: Section 1 - Age Data Integration
4. **Cell 69-70**: Section 2 - scVI Batch Correction
5. **Cell 71-72**: Section 3 - Bayesian GPLVM
6. **Cell 73-74**: Section 4 - AMH Calibration
7. **Cell 75-76**: Section 5 - Risk Stratification
8. **Cell 77-78**: Section 6 - Cross-Study Validation
9. **Cell 79-80**: Section 7 - Final Results Integration

## Common Issues and Solutions

### Issue: "ModuleNotFoundError: No module named 'scvi'"
**Solution**: 
```python
!pip install scvi-tools
```

### Issue: "ModuleNotFoundError: No module named 'gpflow'"
**Solution**:
```python
!pip install gpflow tensorflow
```

### Issue: "NameError: name 'adata' is not defined"
**Solution**: Run the data loading cells first (Cells 1-64)

### Issue: scVI training fails
**Solution**: The code falls back to PCA representation automatically

### Issue: GPLVM training fails
**Solution**: Check tensorflow/gpflow installation and ensure data is numeric

## Testing

Run the test script to check imports:
```bash
python3 test_upgrade_sections.py
```

## Notes

- All sections include comprehensive error handling
- Missing packages will trigger warnings but won't crash the notebook
- Sections will use fallback methods when primary methods fail
- All outputs are saved to files for later inspection

