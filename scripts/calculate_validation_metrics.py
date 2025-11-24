"""
Calculate Validation Metrics for Model Performance

This script implements:
1. Train/test split performance evaluation
2. Cross-validation (K-fold and Leave-One-Out)
3. Classification metrics (AUC-ROC, sensitivity, specificity)
4. Prediction error metrics (MAE, R²)

Author: Agna Chan, Aniqa Nayim, Rimjhim Singh
Date: November 2025
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, KFold, LeaveOneOut
from sklearn.metrics import roc_auc_score, roc_curve, confusion_matrix
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import pearsonr, spearmanr
import json
import os
import warnings
warnings.filterwarnings('ignore')

def calculate_classification_metrics(y_true, y_pred_proba, threshold=0.5):
    """
    Calculate classification metrics from predicted probabilities.
    
    Parameters
    ----------
    y_true : array-like
        True binary labels
    y_pred_proba : array-like
        Predicted probabilities
    threshold : float
        Classification threshold
    
    Returns
    -------
    metrics : dict
        Dictionary of classification metrics
    """
    y_pred = (y_pred_proba >= threshold).astype(int)
    
    # AUC-ROC
    try:
        auc = roc_auc_score(y_true, y_pred_proba)
        fpr, tpr, thresholds_roc = roc_curve(y_true, y_pred_proba)
    except:
        auc = np.nan
        fpr, tpr, thresholds_roc = np.array([]), np.array([]), np.array([])
    
    # Confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    if cm.size == 4:
        tn, fp, fn, tp = cm.ravel()
        sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    else:
        sensitivity = specificity = ppv = npv = np.nan
    
    return {
        'auc_roc': float(auc) if not np.isnan(auc) else np.nan,
        'sensitivity': float(sensitivity) if not np.isnan(sensitivity) else np.nan,
        'specificity': float(specificity) if not np.isnan(specificity) else np.nan,
        'ppv': float(ppv) if not np.isnan(ppv) else np.nan,
        'npv': float(npv) if not np.isnan(npv) else np.nan,
        'fpr': fpr.tolist() if len(fpr) > 0 else [],
        'tpr': tpr.tolist() if len(tpr) > 0 else [],
        'thresholds': thresholds_roc.tolist() if len(thresholds_roc) > 0 else []
    }

def train_test_validation(df, target_col, feature_col, test_size=0.3, random_state=42):
    """
    Perform train/test split validation.
    
    Parameters
    ----------
    df : DataFrame
        Input data
    target_col : str
        Target variable column name
    feature_col : str
        Feature column name for prediction
    test_size : float
        Proportion of data for testing
    random_state : int
        Random seed
    
    Returns
    -------
    results : dict
        Validation results
    """
    # Filter valid data
    mask = df[target_col].notna() & df[feature_col].notna()
    df_valid = df[mask].copy()
    
    if len(df_valid) < 4:
        return {'error': 'Insufficient data for train/test split'}
    
    X = df_valid[[feature_col]].values
    y = df_valid[target_col].values
    
    # Train/test split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state
    )
    
    # Simple linear model for prediction
    from sklearn.linear_model import LinearRegression
    model = LinearRegression()
    model.fit(X_train, y_train)
    
    # Predictions
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)
    
    # Metrics
    train_corr, train_p = pearsonr(y_train, y_train_pred)
    test_corr, test_p = pearsonr(y_test, y_test_pred)
    
    train_mae = mean_absolute_error(y_train, y_train_pred)
    test_mae = mean_absolute_error(y_test, y_test_pred)
    
    train_r2 = r2_score(y_train, y_train_pred)
    test_r2 = r2_score(y_test, y_test_pred)
    
    return {
        'n_train': len(y_train),
        'n_test': len(y_test),
        'train_correlation': float(train_corr),
        'train_pvalue': float(train_p),
        'train_mae': float(train_mae),
        'train_r2': float(train_r2),
        'test_correlation': float(test_corr),
        'test_pvalue': float(test_p),
        'test_mae': float(test_mae),
        'test_r2': float(test_r2)
    }

def cross_validation_metrics(df, target_col, feature_col, cv_type='kfold', n_splits=5, random_state=42):
    """
    Perform cross-validation.
    
    Parameters
    ----------
    df : DataFrame
        Input data
    target_col : str
        Target variable column name
    feature_col : str
        Feature column name
    cv_type : str
        'kfold' or 'loo' (Leave-One-Out)
    n_splits : int
        Number of folds (for KFold)
    random_state : int
        Random seed
    
    Returns
    -------
    results : dict
        CV results
    """
    # Filter valid data
    mask = df[target_col].notna() & df[feature_col].notna()
    df_valid = df[mask].copy()
    
    if len(df_valid) < 3:
        return {'error': 'Insufficient data for cross-validation'}
    
    X = df_valid[[feature_col]].values
    y = df_valid[target_col].values
    
    # Choose CV method
    if cv_type == 'loo' or len(df_valid) < 10:
        cv = LeaveOneOut()
        n_splits = len(df_valid)
    else:
        cv = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    
    cv_correlations = []
    cv_maes = []
    cv_r2s = []
    
    from sklearn.linear_model import LinearRegression
    
    for train_idx, test_idx in cv.split(X):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]
        
        if len(X_train) < 2 or len(X_test) < 1:
            continue
        
        model = LinearRegression()
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        
        if len(y_test) == 1:
            # For LOO, calculate correlation differently
            corr = np.nan
        else:
            corr, _ = pearsonr(y_test, y_pred)
        
        mae = mean_absolute_error(y_test, y_pred)
        r2 = r2_score(y_test, y_pred)
        
        if not np.isnan(corr):
            cv_correlations.append(corr)
        cv_maes.append(mae)
        cv_r2s.append(r2)
    
    if len(cv_correlations) == 0:
        cv_correlations = [np.nan]
    
    return {
        'cv_type': cv_type,
        'n_splits': n_splits,
        'n_folds_completed': len(cv_correlations),
        'mean_correlation': float(np.nanmean(cv_correlations)),
        'std_correlation': float(np.nanstd(cv_correlations)),
        'mean_mae': float(np.mean(cv_maes)),
        'std_mae': float(np.std(cv_maes)),
        'mean_r2': float(np.mean(cv_r2s)),
        'std_r2': float(np.std(cv_r2s)),
        'cv_consistency': float(np.nanstd(cv_correlations) / np.nanmean(cv_correlations)) if np.nanmean(cv_correlations) != 0 else np.nan
    }

def stage_classification_metrics(df):
    """
    Calculate classification metrics for GV vs MI stage prediction.
    
    Parameters
    ----------
    df : DataFrame
        Input data with 'stage' and prediction columns
    
    Returns
    -------
    results : dict
        Classification metrics
    """
    results = {}
    
    # Filter to GV and MI only
    stage_mask = df['stage'].isin(['GV', 'MI'])
    df_stage = df[stage_mask].copy()
    
    if len(df_stage) < 4:
        return {'error': 'Insufficient data for classification'}
    
    # Binary labels: GV = 1, MI = 0
    df_stage['is_GV'] = (df_stage['stage'] == 'GV').astype(int)
    
    # Try different features for classification
    features_to_test = ['health_score', 'cellular_age_z', 'risk_score']
    
    for feature in features_to_test:
        if feature not in df_stage.columns:
            continue
        
        mask = df_stage[feature].notna() & df_stage['is_GV'].notna()
        df_feat = df_stage[mask].copy()
        
        if len(df_feat) < 4:
            continue
        
        y_true = df_feat['is_GV'].values
        
        # Use feature as probability (normalize to 0-1)
        y_pred_proba = df_feat[feature].values
        y_pred_proba = (y_pred_proba - y_pred_proba.min()) / (y_pred_proba.max() - y_pred_proba.min() + 1e-8)
        
        # For health_score, higher = GV (invert for cellular_age_z)
        if feature == 'cellular_age_z' or feature == 'risk_score':
            y_pred_proba = 1 - y_pred_proba
        
        metrics = calculate_classification_metrics(y_true, y_pred_proba)
        results[feature] = metrics
    
    return results

def main():
    """Main function to calculate validation metrics."""
    
    # Load data
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
        print("Error: clinical_decision_framework_final.csv not found")
        return
    
    df = pd.read_csv(clinical_csv, index_col=0)
    print(f"Loaded clinical data: {len(df)} samples")
    
    results = {}
    
    # 1. Train/Test Split Validation
    print("\n1. Train/Test Split Validation")
    
    # Health score prediction from cellular age
    if 'health_score' in df.columns and 'cellular_age_z' in df.columns:
        print("\n1.1 Health Score Prediction (from Cellular Age Z)")
        tt_results = train_test_validation(df, 'health_score', 'cellular_age_z')
        if 'error' not in tt_results:
            results['train_test_health_score'] = tt_results
            print(f"  Test set correlation: r = {tt_results['test_correlation']:.3f}, p = {tt_results['test_pvalue']:.4f}")
            print(f"  Test set MAE: {tt_results['test_mae']:.3f}")
            print(f"  Test set R²: {tt_results['test_r2']:.3f}")
    
    # Age prediction from cellular age
    if 'age' in df.columns and 'cellular_age_z' in df.columns:
        print("\n1.2 Age Prediction (from Cellular Age Z)")
        tt_results = train_test_validation(df, 'age', 'cellular_age_z')
        if 'error' not in tt_results:
            results['train_test_age'] = tt_results
            print(f"  Test set correlation: r = {tt_results['test_correlation']:.3f}, p = {tt_results['test_pvalue']:.4f}")
            print(f"  Test set MAE: {tt_results['test_mae']:.2f} years")
            print(f"  Test set R²: {tt_results['test_r2']:.3f}")
    
    # 2. Cross-Validation
    print("\n2. Cross-Validation")
    
    # Health score CV
    if 'health_score' in df.columns and 'cellular_age_z' in df.columns:
        print("\n2.1 Health Score Cross-Validation")
        cv_type = 'loo' if len(df) < 10 else 'kfold'
        cv_results = cross_validation_metrics(df, 'health_score', 'cellular_age_z', cv_type=cv_type)
        if 'error' not in cv_results:
            results['cv_health_score'] = cv_results
            print(f"  CV type: {cv_results['cv_type']}")
            print(f"  CV correlation: r = {cv_results['mean_correlation']:.3f} ± {cv_results['std_correlation']:.3f}")
            print(f"  CV MAE: {cv_results['mean_mae']:.3f} ± {cv_results['std_mae']:.3f}")
            print(f"  CV R²: {cv_results['mean_r2']:.3f} ± {cv_results['std_r2']:.3f}")
            if not np.isnan(cv_results['cv_consistency']):
                print(f"  CV consistency (CV): {cv_results['cv_consistency']:.3f}")
    
    # Age CV
    if 'age' in df.columns and 'cellular_age_z' in df.columns:
        print("\n2.2 Age Cross-Validation")
        cv_type = 'loo' if len(df) < 10 else 'kfold'
        cv_results = cross_validation_metrics(df, 'age', 'cellular_age_z', cv_type=cv_type)
        if 'error' not in cv_results:
            results['cv_age'] = cv_results
            print(f"  CV type: {cv_results['cv_type']}")
            print(f"  CV correlation: r = {cv_results['mean_correlation']:.3f} ± {cv_results['std_correlation']:.3f}")
            print(f"  CV MAE: {cv_results['mean_mae']:.2f} ± {cv_results['std_mae']:.2f} years")
            print(f"  CV R²: {cv_results['mean_r2']:.3f} ± {cv_results['std_r2']:.3f}")
    
    # 3. Classification Metrics (GV vs MI)
    print("\n3. Stage Classification Metrics (GV vs MI)")
    
    # Load sample metadata if available
    sample_csv = os.path.join(os.path.dirname(clinical_csv), 'sample_metadata_with_age.csv')
    if os.path.exists(sample_csv):
        sample_df = pd.read_csv(sample_csv)
        if 'sample' in sample_df.columns and 'stage' in sample_df.columns:
            df = df.merge(sample_df[['sample', 'stage']], 
                         left_index=True, right_on='sample', how='left')
            df.set_index('sample', inplace=True, drop=False)
    
    if 'stage' in df.columns:
        class_results = stage_classification_metrics(df)
        if 'error' not in class_results:
            results['classification'] = class_results
            
            for feature, metrics in class_results.items():
                print(f"\n3.{list(class_results.keys()).index(feature) + 1} Using {feature}:")
                if not np.isnan(metrics['auc_roc']):
                    print(f"  AUC-ROC: {metrics['auc_roc']:.3f}")
                print(f"  Sensitivity: {metrics['sensitivity']:.3f}")
                print(f"  Specificity: {metrics['specificity']:.3f}")
                print(f"  PPV: {metrics['ppv']:.3f}")
                print(f"  NPV: {metrics['npv']:.3f}")
    
    # Save results
    output_dir = os.path.dirname(clinical_csv)
    output_file = os.path.join(output_dir, 'validation_metrics.json')
    
    # Convert numpy types to native Python types
    def convert_to_native(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj) if not np.isnan(obj) else None
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
    print("\nValidation metrics calculation complete.")

if __name__ == '__main__':
    main()

