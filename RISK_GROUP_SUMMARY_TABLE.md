# Table 1. Oocyte Risk Group Summary (Model Outputs)

This table summarizes **EVERYTHING your model outputs**:
- cellular age
- uncertainty
- health score
- clinical risk group

**PERFECT for ADS** because it is interpretable, biological, and summarizes your findings without repeating performance metrics.

| Risk Group | n | Cellular Age (z) | Health Score | Uncertainty (σ) | % of Cells | Biological Interpretation |
|------------|---|------------------|--------------|-----------------|------------|---------------------------|
| **Low Risk** | 13 | Young (0.0-0.3) | 70-75 | Medium (1181 ± 681) | 65.0% | Mostly healthy GV oocytes, strong viability |
| **Moderate Risk** | 6 | Mid (0.3-0.6) | 66-69 | Medium (1969 ± 344) | 30.0% | Transition from GV to MI; early aging signals |
| **High Risk** | 1 | Older (0.6-1.0) | 40-40 | High (15.9K, σ×10³) | 5.0% | Aging MI oocytes; urgent intervention window |

**Key Interpretation**: The risk stratification reveals three distinct oocyte quality trajectories, with 65% showing resilient aging patterns (low risk), 30% in transition (moderate risk), and 5% requiring urgent intervention (high risk). Higher uncertainty values (scaled by 10³) in high-risk cells reflect increased biological variability and decreased prediction confidence, consistent with accelerated aging phenotypes.



