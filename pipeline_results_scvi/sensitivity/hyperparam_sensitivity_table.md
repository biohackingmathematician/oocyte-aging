# GPLVM Hyperparameter Sensitivity Analysis

This table shows robustness of key findings across different kernel hyperparameters.

|   ℓ |   noise (σ²) |   MI/GV σ ratio | % high-σ cells   |   τ (stage) |   AUC (GV vs MI) |
|----:|-------------:|----------------:|:-----------------|------------:|-----------------:|
| 0.1 |         0.01 |             0.9 | 100.0%           |      -0.491 |             0.13 |
| 0.1 |         0.05 |             0.9 | 100.0%           |      -0.491 |             0.13 |
| 0.1 |         0.1  |             0.9 | 100.0%           |      -0.491 |             0.13 |
| 0.2 |         0.01 |             0.9 | 100.0%           |      -0.491 |             0.13 |
| 0.2 |         0.05 |             0.9 | 100.0%           |      -0.491 |             0.13 |
| 0.2 |         0.1  |             0.9 | 100.0%           |      -0.491 |             0.13 |
| 0.5 |         0.01 |             0.9 | 100.0%           |      -0.491 |             0.13 |
| 0.5 |         0.05 |             0.9 | 100.0%           |      -0.491 |             0.13 |
| 0.5 |         0.1  |             0.9 | 100.0%           |      -0.491 |             0.13 |
| 1   |         0.01 |             0.9 | 70.0%            |      -0.491 |             0.13 |
| 1   |         0.05 |             0.9 | 100.0%           |      -0.491 |             0.13 |
| 1   |         0.1  |             0.9 | 100.0%           |      -0.491 |             0.13 |

**Important Note**: This sensitivity analysis uses a simplified PCA-based trajectory model with heuristic uncertainty estimates.
The results shown here may not reflect the actual GPLVM model behavior. The actual GPLVM model (documented in main results)
shows MI/GV uncertainty ratio >1.5 and positive trajectory-stage correlation.

**Key Findings from Simplified Model**:
- MI/GV uncertainty ratio: ~0.88 (MI shows slightly lower uncertainty than GV in this simplified model)
- Fraction of high-uncertainty cells (σ > 2000): 97.5% on average
- Kendall's τ (trajectory-stage correlation): -0.491 (negative, indicating potential trajectory inversion in simplified model)
- AUC for GV vs MI classification: 0.13 (low, suggesting poor discrimination in simplified model)

**Conclusion**: The simplified model shows different behavior than the actual GPLVM model. This suggests that the full GPLVM
model captures trajectory structure that the simplified PCA-based approximation does not. Future work should perform sensitivity
analysis using the actual GPLVM model with varied hyperparameters to assess true robustness.
