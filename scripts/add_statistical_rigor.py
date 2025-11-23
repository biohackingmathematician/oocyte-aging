"""
Add Statistical Rigor Metrics

This script implements:
1. Multiple testing correction (FDR)
2. Bootstrap confidence intervals
3. Effect size calculations (Cohen's d)

Author: Agna Chan, Aniqa Nayim, Rimjhim Singh
Date: November 2025
"""

import pandas as pd
import numpy as np
from scipy.stats import bootstrap
from statsmodels.stats.multitest import multipletests
import json
import os

def calculate_confidence_interval(data, statistic_func, n_resamples=1000, confidence_level=0.95):
    """
    Calculate bootstrap confidence interval for a statistic using manual bootstrap.
    
    Parameters
    ----------
    data : tuple or array-like
        Input data (tuple for multiple arrays)
    statistic_func : callable
        Function to compute statistic (takes data tuple as argument)
    n_resamples : int
        Number of bootstrap resamples
    confidence_level : float
        Confidence level (default 0.95)
    
    Returns
    -------
    ci : tuple
        (lower_bound, upper_bound)
    """
    try:
        np.random.seed(42)
        
        # Ensure data is a tuple
        if isinstance(data, tuple):
            data_tuple = data
        else:
            data_tuple = (data,)
        
        n = len(data_tuple[0])
        bootstrap_stats = []
        
        for _ in range(n_resamples):
            # Resample with replacement
            indices = np.random.choice(n, size=n, replace=True)
            resampled_data = tuple(arr[indices] for arr in data_tuple)
            stat = statistic_func(resampled_data)
            if not np.isnan(stat):
                bootstrap_stats.append(stat)
        
        if len(bootstrap_stats) == 0:
            return (np.nan, np.nan)
        
        alpha = 1 - confidence_level
        lower_percentile = (alpha / 2) * 100
        upper_percentile = (1 - alpha / 2) * 100
        
        ci_lower = np.percentile(bootstrap_stats, lower_percentile)
        ci_upper = np.percentile(bootstrap_stats, upper_percentile)
        
        return (float(ci_lower), float(ci_upper))
    except Exception as e:
        print(f"Warning: Bootstrap failed: {e}")
        return (np.nan, np.nan)

def correlation_with_ci(x, y, n_resamples=1000):
    """
    Calculate Pearson correlation with bootstrap confidence interval.
    
    Parameters
    ----------
    x, y : array-like
        Input arrays
    n_resamples : int
        Number of bootstrap resamples
    
    Returns
    -------
    r : float
        Correlation coefficient
    ci : tuple
        (lower_bound, upper_bound)
    """
    def corr_statistic(data):
        # data is a tuple of (x_resampled, y_resampled)
        x_data, y_data = data
        if len(x_data) < 2 or len(y_data) < 2:
            return np.nan
        corr = np.corrcoef(x_data, y_data)
        if corr.size == 1:
            return np.nan
        return float(corr[0, 1])
    
    r = float(np.corrcoef(x, y)[0, 1])
    ci = calculate_confidence_interval((x, y), corr_statistic, n_resamples)
    return r, ci

def apply_fdr_correction(p_values, alpha=0.05, method='fdr_bh'):
    """
    Apply FDR correction to p-values.
    
    Parameters
    ----------
    p_values : array-like
        Array of p-values
    alpha : float
        Significance level
    method : str
        Correction method ('fdr_bh' for Benjamini-Hochberg)
    
    Returns
    -------
    rejected : array
        Boolean array of rejected hypotheses
    p_adjusted : array
        FDR-adjusted p-values
    """
    p_array = np.array(p_values)
    valid_mask = ~np.isnan(p_array)
    
    if valid_mask.sum() == 0:
        return np.array([]), np.array([])
    
    rejected = np.zeros(len(p_array), dtype=bool)
    p_adjusted = np.full(len(p_array), np.nan)
    
    if valid_mask.sum() > 0:
        _, p_adj, _, _ = multipletests(p_array[valid_mask], alpha=alpha, method=method)
        rejected[valid_mask] = p_adj < alpha
        p_adjusted[valid_mask] = p_adj
    
    return rejected, p_adjusted

def calculate_cohens_d(group1, group2):
    """
    Calculate Cohen's d effect size.
    
    Parameters
    ----------
    group1, group2 : array-like
        Two groups to compare
    
    Returns
    -------
    d : float
        Cohen's d effect size
    interpretation : str
        Effect size interpretation
    """
    g1 = np.array(group1)
    g2 = np.array(group2)
    
    g1 = g1[~np.isnan(g1)]
    g2 = g2[~np.isnan(g2)]
    
    if len(g1) == 0 or len(g2) == 0:
        return np.nan, "Insufficient data"
    
    mean1, mean2 = np.mean(g1), np.mean(g2)
    std1, std2 = np.std(g1, ddof=1), np.std(g2, ddof=1)
    
    # Pooled standard deviation
    n1, n2 = len(g1), len(g2)
    pooled_std = np.sqrt(((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2))
    
    if pooled_std == 0:
        return np.nan, "No variance"
    
    d = (mean1 - mean2) / pooled_std
    
    # Interpretation
    abs_d = abs(d)
    if abs_d < 0.2:
        interpretation = "Small effect"
    elif abs_d < 0.5:
        interpretation = "Medium effect"
    elif abs_d < 0.8:
        interpretation = "Large effect"
    else:
        interpretation = "Very large effect"
    
    return d, interpretation

def main():
    """Main function to add statistical rigor metrics."""
    
    # Load data
    # Try multiple possible paths
    possible_paths = [
        '../data/clinical_decision_framework_final.csv',
        'data/clinical_decision_framework_final.csv',
        os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'clinical_decision_framework_final.csv')
    ]
    
    clinical_csv = None
    for path in possible_paths:
        if os.path.exists(path):
            clinical_csv = path
            break
    
    if clinical_csv is None:
        print("Error: clinical_decision_framework_final.csv not found in any expected location")
        print(f"Tried: {possible_paths}")
        return
    
    df = pd.read_csv(clinical_csv, index_col=0)
    print(f"Loaded clinical data: {len(df)} samples")
    
    results = {}
    
    # 1. Multiple Testing Correction
    print("\n1. Multiple Testing Correction")
    
    # Collect all p-values from correlations
    p_values = []
    correlation_names = []
    
    if 'cellular_age_z' in df.columns:
        if 'age' in df.columns:
            age_mask = df['age'].notna() & df['cellular_age_z'].notna()
            if age_mask.sum() >= 3:
                from scipy.stats import pearsonr
                _, pval = pearsonr(df.loc[age_mask, 'age'], 
                                  df.loc[age_mask, 'cellular_age_z'])
                p_values.append(pval)
                correlation_names.append('z_vs_age')
        
        if 'health_score' in df.columns:
            hs_mask = df['health_score'].notna() & df['cellular_age_z'].notna()
            if hs_mask.sum() >= 3:
                from scipy.stats import pearsonr
                _, pval = pearsonr(df.loc[hs_mask, 'cellular_age_z'],
                                  df.loc[hs_mask, 'health_score'])
                p_values.append(pval)
                correlation_names.append('z_vs_health_score')
    
    if len(p_values) > 0:
        rejected, p_adjusted = apply_fdr_correction(p_values)
        results['multiple_testing'] = {
            'n_tests': len(p_values),
            'n_significant_uncorrected': (np.array(p_values) < 0.05).sum(),
            'n_significant_fdr': rejected.sum(),
            'fdr_adjusted_pvalues': {name: float(p_adj) for name, p_adj in 
                                     zip(correlation_names, p_adjusted) if not np.isnan(p_adj)}
        }
        print(f"Tests performed: {len(p_values)}")
        print(f"Significant (uncorrected): {(np.array(p_values) < 0.05).sum()}")
        print(f"Significant (FDR < 0.05): {rejected.sum()}")
        for name, p_adj in zip(correlation_names, p_adjusted):
            if not np.isnan(p_adj):
                print(f"  {name}: p_adj = {p_adj:.4f}")
    
    # 2. Confidence Intervals
    print("\n2. Confidence Intervals")
    
    ci_results = {}
    
    if 'cellular_age_z' in df.columns:
        if 'age' in df.columns:
            age_mask = df['age'].notna() & df['cellular_age_z'].notna()
            if age_mask.sum() >= 3:
                x = df.loc[age_mask, 'age'].values
                y = df.loc[age_mask, 'cellular_age_z'].values
                r, ci = correlation_with_ci(x, y)
                ci_results['z_vs_age_correlation'] = {
                    'r': float(r),
                    'ci_lower': float(ci[0]),
                    'ci_upper': float(ci[1])
                }
                print(f"Z vs Age correlation: r = {r:.3f} (95% CI: {ci[0]:.3f}-{ci[1]:.3f})")
        
        if 'health_score' in df.columns:
            hs_mask = df['health_score'].notna() & df['cellular_age_z'].notna()
            if hs_mask.sum() >= 3:
                x = df.loc[hs_mask, 'cellular_age_z'].values
                y = df.loc[hs_mask, 'health_score'].values
                r, ci = correlation_with_ci(x, y)
                ci_results['z_vs_health_score_correlation'] = {
                    'r': float(r),
                    'ci_lower': float(ci[0]),
                    'ci_upper': float(ci[1])
                }
                print(f"Z vs Health Score correlation: r = {r:.3f} (95% CI: {ci[0]:.3f}-{ci[1]:.3f})")
    
    # Health score means by stage
    if 'health_score' in df.columns and 'stage' in df.columns:
        stage_mask = df['stage'].isin(['GV', 'MI'])
        if stage_mask.sum() > 0:
            gv_scores = df[df['stage'] == 'GV']['health_score'].dropna().values
            mi_scores = df[df['stage'] == 'MI']['health_score'].dropna().values
            
            if len(gv_scores) > 0:
                def mean_statistic(data):
                    return float(np.mean(data[0]))
                ci_gv = calculate_confidence_interval((gv_scores,), mean_statistic)
                ci_results['health_score_gv'] = {
                    'mean': float(np.mean(gv_scores)),
                    'ci_lower': float(ci_gv[0]),
                    'ci_upper': float(ci_gv[1])
                }
                print(f"GV health score: {np.mean(gv_scores):.1f} (95% CI: {ci_gv[0]:.1f}-{ci_gv[1]:.1f})")
            
            if len(mi_scores) > 0:
                def mean_statistic(data):
                    return float(np.mean(data[0]))
                ci_mi = calculate_confidence_interval((mi_scores,), mean_statistic)
                ci_results['health_score_mi'] = {
                    'mean': float(np.mean(mi_scores)),
                    'ci_lower': float(ci_mi[0]),
                    'ci_upper': float(ci_mi[1])
                }
                print(f"MI health score: {np.mean(mi_scores):.1f} (95% CI: {ci_mi[0]:.1f}-{ci_mi[1]:.1f})")
    
    results['confidence_intervals'] = ci_results
    
    # 3. Effect Sizes
    print("\n3. Effect Sizes (Cohen's d)")
    
    effect_sizes = {}
    
    if 'health_score' in df.columns and 'stage' in df.columns:
        stage_mask = df['stage'].isin(['GV', 'MI'])
        if stage_mask.sum() > 0:
            gv_scores = df[df['stage'] == 'GV']['health_score'].dropna().values
            mi_scores = df[df['stage'] == 'MI']['health_score'].dropna().values
            
            if len(gv_scores) > 0 and len(mi_scores) > 0:
                d, interpretation = calculate_cohens_d(gv_scores, mi_scores)
                effect_sizes['gv_vs_mi_health_score'] = {
                    'cohens_d': float(d),
                    'interpretation': interpretation
                }
                print(f"GV vs MI health score: Cohen's d = {d:.3f} ({interpretation})")
    
    results['effect_sizes'] = effect_sizes
    
    # Save results
    output_dir = os.path.dirname(clinical_csv)
    output_file = os.path.join(output_dir, 'statistical_rigor_metrics.json')
    
    # Convert numpy types to native Python types for JSON serialization
    def convert_to_native(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {key: convert_to_native(value) for key, value in obj.items()}
        elif isinstance(obj, list):
            return [convert_to_native(item) for item in obj]
        return obj
    
    results_native = convert_to_native(results)
    with open(output_file, 'w') as f:
        json.dump(results_native, f, indent=2)
    
    print(f"\nResults saved to {output_file}")
    print("\nStatistical rigor metrics calculation complete.")

if __name__ == '__main__':
    main()

