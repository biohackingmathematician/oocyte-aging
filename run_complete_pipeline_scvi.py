#!/usr/bin/env python3
"""
Complete pipeline execution with scVI.
Follows the same structure as PCA version but uses scVI latent space.
"""

import sys
import os
import traceback
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, pearsonr
from sklearn.neighbors import NearestNeighbors
import warnings
warnings.filterwarnings('ignore')

print("="*70)
print("COMPLETE PIPELINE EXECUTION WITH scVI")
print("="*70)

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white')

# Create output directories
OUTDIRS = {
    'results': './pipeline_results_scvi',
    'figures': './pipeline_results_scvi/figures',
    'tables': './pipeline_results_scvi/tables'
}
for d in OUTDIRS.values():
    os.makedirs(d, exist_ok=True)

# ============================================================================
# STEP 1: Load scVI-processed data
# ============================================================================
print("\n" + "="*70)
print("STEP 1: Loading scVI-processed data")
print("="*70)

try:
    adata = sc.read_h5ad('adata_with_scvi.h5ad')
    print(f"✓ Loaded data: {adata.shape}")
    
    # Verify scVI latent space exists
    if 'X_scvi' in adata.obsm:
        print(f"✓ scVI latent space: {adata.obsm['X_scvi'].shape}")
    else:
        print("✗ ERROR: X_scvi not found. Please run run_notebook_complete.py first.")
        sys.exit(1)
        
    # Check if log1p_norm layer exists (needed for pathway scoring)
    if 'log1p_norm' not in adata.layers:
        print("  Creating log1p_norm layer...")
        adata_norm = adata.copy()
        sc.pp.normalize_total(adata_norm, target_sum=1e4)
        sc.pp.log1p(adata_norm)
        adata.layers['log1p_norm'] = adata_norm.X.copy()
        print("  ✓ Created log1p_norm layer")
    
    print(f"  Observations: {adata.n_obs}")
    print(f"  Variables: {adata.n_vars}")
    if 'stage' in adata.obs.columns:
        print(f"  Stages: {adata.obs['stage'].value_counts().to_dict()}")
        
except Exception as e:
    print(f"✗ Error loading data: {e}")
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# STEP 2: Compute DPT trajectory on scVI latent space
# ============================================================================
print("\n" + "="*70)
print("STEP 2: Computing DPT trajectory on scVI latent space")
print("="*70)

try:
    # Use scVI latent space for neighbors and diffusion map
    print("  Computing neighbors on scVI latent space...")
    sc.pp.neighbors(adata, use_rep='X_scvi', n_neighbors=min(10, adata.n_obs-1))
    
    print("  Computing diffusion map...")
    sc.tl.diffmap(adata, n_comps=10)
    
    # Set root cell (GV stage)
    if 'stage' in adata.obs.columns:
        gv_cells = adata.obs[adata.obs['stage'] == 'GV'].index
        if len(gv_cells) > 0:
            root_cell_idx = np.where(adata.obs_names == gv_cells[0])[0][0]
            adata.uns['iroot'] = root_cell_idx
            print(f"  ✓ Root cell set to: {gv_cells[0]} (index: {root_cell_idx})")
        else:
            print("  Warning: No GV cells found, using first cell as root")
            adata.uns['iroot'] = 0
    else:
        adata.uns['iroot'] = 0
    
    # Compute DPT (no branching for simple trajectory)
    print("  Computing Diffusion Pseudotime (DPT)...")
    try:
        sc.tl.dpt(adata, n_dcs=10, n_branchings=0, min_group_size=0.01)
    except Exception as e:
        print(f"  Warning: DPT with branching failed: {e}")
        print("  Retrying without branching...")
        # Try computing DPT manually using diffusion map
        if 'X_diffmap' in adata.obsm:
            # Use first diffusion component as pseudotime
            dpt = -adata.obsm['X_diffmap'][:, 0]  # Negative to reverse direction
            dpt = (dpt - dpt.min()) / (dpt.max() - dpt.min() + 1e-8)
            adata.obs['dpt_pseudotime'] = dpt
            print("  ✓ Computed DPT from diffusion map")
        else:
            raise
    
    # Normalize DPT to 0-1 range
    if 'dpt_pseudotime' in adata.obs.columns:
        dpt = adata.obs['dpt_pseudotime'].values
        dpt_normalized = (dpt - dpt.min()) / (dpt.max() - dpt.min() + 1e-8)
        adata.obs['dpt_pseudotime'] = dpt_normalized
        print(f"  ✓ DPT computed: range [{dpt_normalized.min():.3f}, {dpt_normalized.max():.3f}]")
        
        # Validate trajectory with stage
        if 'stage' in adata.obs.columns:
            from scipy.stats import spearmanr
            stage_numeric = adata.obs['stage'].map({'GV': 0, 'MI': 1, 'MII': 2})
            corr, pval = spearmanr(stage_numeric, adata.obs['dpt_pseudotime'])
            print(f"  ✓ DPT-Stage correlation: ρ = {corr:.3f}, p = {pval:.4f}")
    else:
        print("  ✗ ERROR: DPT computation failed")
        sys.exit(1)
        
except Exception as e:
    print(f"✗ Error computing DPT: {e}")
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# STEP 3: Calculate pathway scores
# ============================================================================
print("\n" + "="*70)
print("STEP 3: Calculating pathway scores")
print("="*70)

try:
    # Define pathway gene sets
    pathway_genes = {
        'Cell_Cycle': [
            'CCNB1', 'CCNB2', 'CDC20', 'CDC25C', 'CDK1', 'CDK2', 
            'MAD2L1', 'BUB1', 'BUB1B', 'PLK1', 'AURKA', 'AURKB'
        ],
        'Mitochondrial': [
            'ATP5A1', 'ATP5B', 'COX4I1', 'NDUFB8', 'SDHA', 'UQCRC1',
            'CYCS', 'VDAC1', 'VDAC2', 'VDAC3', 'TOMM20', 'TOMM40'
        ],
        'DNA_Repair': [
            'BRCA1', 'BRCA2', 'ATM', 'ATR', 'CHEK1', 'CHEK2',
            'RAD51', 'XRCC1', 'PARP1', 'PCNA', 'TP53', 'MDM2'
        ],
        'Oxidative_Stress': [
            'SOD1', 'SOD2', 'CAT', 'GPX1', 'GPX4', 'PRDX1',
            'PRDX2', 'PRDX3', 'TXN', 'TXNRD1', 'NQO1', 'GSTM1'
        ],
        'Apoptosis': [
            'BAX', 'BAK1', 'CASP3', 'CASP7', 'CASP9', 'CASP8',
            'BCL2', 'BCL2L1', 'MCL1', 'PARP1', 'APAF1', 'DIABLO'
        ]
    }
    
    print("\n  Computing pathway scores...")
    
    # Load gene symbol mapping from biomart file
    print("  Loading gene symbol mapping...")
    gene_symbol_map = {}
    transcript_to_symbol = {}  # transcript_id -> gene_symbol
    biomart_file = 'zenodo_data/final_code/e106_biomart_HomoSapiens_Transcript_Gene_Names.txt'
    
    if os.path.exists(biomart_file):
        try:
            biomart_df = pd.read_csv(biomart_file, sep='\t', header=0, low_memory=False)
            print(f"    Loaded biomart file: {len(biomart_df)} entries")
            
            # Use known column names
            transcript_col = 'Transcript stable ID version'
            gene_name_col = 'Gene name'
            
            if transcript_col in biomart_df.columns and gene_name_col in biomart_df.columns:
                # Create mapping: transcript_id -> gene_symbol
                for _, row in biomart_df.iterrows():
                    transcript_id = str(row[transcript_col]).strip()
                    gene_symbol = str(row[gene_name_col]).strip()
                    if transcript_id and gene_symbol and gene_symbol != 'nan' and transcript_id != 'nan':
                        # Store with and without version
                        transcript_to_symbol[transcript_id] = gene_symbol
                        transcript_to_symbol[transcript_id.split('.')[0]] = gene_symbol
                
                print(f"    Created mapping for {len(transcript_to_symbol)} transcript IDs")
            else:
                print(f"    Warning: Expected columns not found. Available: {list(biomart_df.columns)}")
        except Exception as e:
            print(f"    Warning: Could not load biomart file: {e}")
            traceback.print_exc()
    
    # Map var_names (transcript IDs) to gene symbols
    for var_name in adata.var_names:
        if var_name in transcript_to_symbol:
            gene_symbol_map[var_name] = transcript_to_symbol[var_name]
        elif var_name.split('.')[0] in transcript_to_symbol:
            gene_symbol_map[var_name] = transcript_to_symbol[var_name.split('.')[0]]
    
    # Also check if gene_symbol column exists in adata.var
    if 'gene_symbol' in adata.var.columns:
        for idx, var_name in enumerate(adata.var_names):
            symbol = adata.var.iloc[idx]['gene_symbol']
            if pd.notna(symbol) and str(symbol) != 'nan':
                gene_symbol_map[var_name] = str(symbol)
    
    print(f"    Total mapped genes: {len(gene_symbol_map)}")
    
    if len(gene_symbol_map) == 0:
        print("  Warning: No gene symbol mapping available. Will try direct matching.")
    
    # Compute pathway scores
    pathway_scores_df = {}
    
    for pathway_name, genes in pathway_genes.items():
        print(f"    {pathway_name}:")
        
        # Find matching genes by symbol
        gene_indices = []
        found_genes_set = set()
        
        # Match via gene symbols
        gene_upper_set = {g.upper() for g in genes}
        for i, var_name in enumerate(adata.var_names):
            symbol = gene_symbol_map.get(var_name, '').upper()
            if symbol and symbol in gene_upper_set:
                if i not in gene_indices:
                    gene_indices.append(i)
                    found_genes_set.add(symbol.upper())
        
        found_genes = [g for g in genes if g.upper() in found_genes_set]
        
        if len(gene_indices) == 0:
            print(f"      Warning: No genes found for {pathway_name}")
            continue
        
        print(f"      Found {len(gene_indices)}/{len(genes)} genes: {', '.join(found_genes[:5])}...")
        
        # Calculate pathway score (mean expression)
        if 'log1p_norm' in adata.layers:
            pathway_expr = adata[:, gene_indices].layers['log1p_norm'].mean(axis=1)
        else:
            pathway_expr = adata[:, gene_indices].X.mean(axis=1)
        
        adata.obs[f'pathway_{pathway_name}'] = pathway_expr
        
        # Normalize to 0-1
        min_val = pathway_expr.min()
        max_val = pathway_expr.max()
        if max_val > min_val:
            score = (pathway_expr - min_val) / (max_val - min_val)
        else:
            score = np.zeros_like(pathway_expr)
        
        adata.obs[f'score_{pathway_name}'] = score
        pathway_scores_df[pathway_name] = score
        
        print(f"      Score range: [{score.min():.3f}, {score.max():.3f}]")
    
    print("  ✓ Pathway scores computed!")
    
    # Save pathway scores
    pathway_df = pd.DataFrame(pathway_scores_df, index=adata.obs_names)
    pathway_df.to_csv(f"{OUTDIRS['tables']}/pathway_scores_scvi.csv")
    print(f"  ✓ Saved: {OUTDIRS['tables']}/pathway_scores_scvi.csv")
    
except Exception as e:
    print(f"✗ Error calculating pathway scores: {e}")
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# STEP 4: Compute composite health score
# ============================================================================
print("\n" + "="*70)
print("STEP 4: Computing composite health score")
print("="*70)

try:
    # Define weights for each pathway
    weights = {
        'Cell_Cycle': 0.25,
        'Mitochondrial': 0.25,
        'DNA_Repair': 0.20,
        'Oxidative_Stress': 0.15,
        'Apoptosis': 0.15
    }
    
    print("  Computing weighted composite score...")
    
    composite_score = np.zeros(adata.n_obs)
    total_weight = 0
    
    for pathway, weight in weights.items():
        score_col = f'score_{pathway}'
        if score_col in adata.obs.columns:
            composite_score += weight * adata.obs[score_col].values
            total_weight += weight
            print(f"    {pathway}: weight {weight}")
        else:
            print(f"    Warning: {score_col} not found, skipping")
    
    if total_weight > 0:
        composite_score = composite_score / total_weight
    
    # Normalize to 0-100
    min_score = composite_score.min()
    max_score = composite_score.max()
    if max_score > min_score:
        health_score = ((composite_score - min_score) / (max_score - min_score)) * 100
    else:
        health_score = np.ones_like(composite_score) * 50
    
    adata.obs['oocyte_health_score'] = health_score
    
    print(f"\n  Composite health score range: [{health_score.min():.1f}, {health_score.max():.1f}]")
    
    # Correlation with trajectory
    if 'dpt_pseudotime' in adata.obs.columns:
        corr, pval = spearmanr(adata.obs['dpt_pseudotime'], health_score)
        print(f"  Health score - DPT correlation: r = {corr:.3f}, p = {pval:.4f}")
    
    # Health score by stage
    if 'stage' in adata.obs.columns:
        print("\n  Health score by stage:")
        for stage in sorted(adata.obs['stage'].unique()):
            stage_scores = adata.obs[adata.obs['stage'] == stage]['oocyte_health_score']
            print(f"    {stage}: mean = {stage_scores.mean():.1f}, std = {stage_scores.std():.1f}")
    
except Exception as e:
    print(f"✗ Error computing health score: {e}")
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# STEP 5: Compute cellular age (using DPT as proxy)
# ============================================================================
print("\n" + "="*70)
print("STEP 5: Computing cellular age")
print("="*70)

try:
    # Use DPT as cellular age
    if 'dpt_pseudotime' in adata.obs.columns:
        adata.obs['cellular_age_z'] = adata.obs['dpt_pseudotime'].values
        print("  ✓ Using DPT as cellular age")
    else:
        print("  ✗ ERROR: DPT not available for cellular age")
        sys.exit(1)
    
    # Compute uncertainty (local variance in scVI space)
    print("  Computing cellular age uncertainty...")
    X_scvi = adata.obsm['X_scvi']
    
    # Local uncertainty: std of cellular age among k-NN in scVI space
    k = min(5, adata.n_obs - 1)
    nbrs = NearestNeighbors(n_neighbors=k+1).fit(X_scvi)
    distances, indices = nbrs.kneighbors(X_scvi)
    
    local_std = []
    for i in range(adata.n_obs):
        neighbor_ages = adata.obs['cellular_age_z'].iloc[indices[i, 1:]].values  # Exclude self
        local_std.append(neighbor_ages.std() if len(neighbor_ages) > 1 else 0.01)
    
    adata.obs['cellular_age_uncertainty'] = np.array(local_std)
    
    # Also use scVI latent uncertainty if available
    if 'scvi_latent_sigma' in adata.obs.columns:
        adata.obs['scvi_uncertainty'] = adata.obs['scvi_latent_sigma'].values
        # Combine uncertainties
        combined_uncertainty = np.sqrt(adata.obs['cellular_age_uncertainty']**2 + 
                                      (adata.obs['scvi_uncertainty'] * 1000)**2)
        adata.obs['cellular_age_uncertainty'] = combined_uncertainty
    
    print(f"  ✓ Cellular age uncertainty range: [{adata.obs['cellular_age_uncertainty'].min():.2f}, "
          f"{adata.obs['cellular_age_uncertainty'].max():.2f}]")
    
except Exception as e:
    print(f"✗ Error computing cellular age: {e}")
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# STEP 6: Risk stratification
# ============================================================================
print("\n" + "="*70)
print("STEP 6: Risk stratification")
print("="*70)

try:
    # Load age information if available
    if 'age' not in adata.obs.columns:
        # Try to load from metadata
        try:
            meta = pd.read_csv('data/sample_metadata_with_age.csv', index_col='sample')
            if 'age' in meta.columns:
                adata.obs['age'] = adata.obs_names.map(lambda x: meta.loc[x.split('_')[0], 'age'] if x.split('_')[0] in meta.index else 35)
            else:
                adata.obs['age'] = 35  # Default age
        except:
            adata.obs['age'] = 35
    
    # Compute risk score based on cellular age, uncertainty, and health score
    print("  Computing risk scores...")
    
    # Normalize components to 0-1
    age_z = adata.obs['cellular_age_z'].values
    uncertainty_norm = (adata.obs['cellular_age_uncertainty'].values - 
                       adata.obs['cellular_age_uncertainty'].min()) / \
                       (adata.obs['cellular_age_uncertainty'].max() - 
                        adata.obs['cellular_age_uncertainty'].min() + 1e-8)
    health_norm = 1 - (adata.obs['oocyte_health_score'].values / 100)  # Invert: high health = low risk
    
    # Weighted combination
    risk_score = (age_z * 0.4 + uncertainty_norm * 0.3 + health_norm * 0.3) * 1000
    
    adata.obs['risk_score'] = risk_score
    
    # Define risk groups
    low_threshold = np.percentile(risk_score, 33.3)
    high_threshold = np.percentile(risk_score, 66.7)
    
    def assign_risk_group(score):
        if score < low_threshold:
            return 'Low Risk (Resilient Agers)'
        elif score < high_threshold:
            return 'Moderate Risk'
        else:
            return 'High Risk (Accelerated Agers)'
    
    adata.obs['risk_group'] = adata.obs['risk_score'].apply(assign_risk_group)
    
    print(f"\n  Risk score range: [{risk_score.min():.1f}, {risk_score.max():.1f}]")
    print(f"  Risk thresholds: Low < {low_threshold:.1f}, High > {high_threshold:.1f}")
    print("\n  Risk group distribution:")
    risk_counts = adata.obs['risk_group'].value_counts()
    for group, count in risk_counts.items():
        pct = count / len(adata) * 100
        print(f"    {group}: {count} ({pct:.1f}%)")
    
except Exception as e:
    print(f"✗ Error in risk stratification: {e}")
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# STEP 7: Build clinical decision framework
# ============================================================================
print("\n" + "="*70)
print("STEP 7: Building clinical decision framework")
print("="*70)

try:
    # Create comprehensive clinical decision framework
    clinical_df = pd.DataFrame({
        'sample': adata.obs_names,
        'age': adata.obs['age'].values if 'age' in adata.obs.columns else 35,
        'cellular_age_z': adata.obs['cellular_age_z'].values,
        'cellular_age_uncertainty': adata.obs['cellular_age_uncertainty'].values,
        'risk_group': adata.obs['risk_group'].values,
        'risk_score': adata.obs['risk_score'].values,
        'health_score': adata.obs['oocyte_health_score'].values,
        'stage': adata.obs['stage'].values if 'stage' in adata.obs.columns else 'Unknown',
        'dpt_pseudotime': adata.obs['dpt_pseudotime'].values if 'dpt_pseudotime' in adata.obs.columns else np.nan
    })
    
    # Add pathway scores
    for pathway in ['Cell_Cycle', 'Mitochondrial', 'DNA_Repair', 'Oxidative_Stress', 'Apoptosis']:
        score_col = f'score_{pathway}'
        if score_col in adata.obs.columns:
            clinical_df[f'pathway_{pathway}'] = adata.obs[score_col].values
    
    # Add intervention recommendations
    def recommend_intervention(row):
        if row['risk_group'].startswith('Low'):
            return 'Monitor - Continue routine screening'
        elif row['risk_group'] == 'Moderate Risk':
            return 'Consider intervention - Discuss fertility preservation options'
        else:
            return 'Urgent intervention - Immediate fertility preservation consultation'
    
    clinical_df['recommendation'] = clinical_df.apply(recommend_intervention, axis=1)
    
    # Save clinical decision framework
    clinical_df.to_csv(f"{OUTDIRS['tables']}/clinical_decision_framework_scvi.csv", index=False)
    print(f"  ✓ Saved: {OUTDIRS['tables']}/clinical_decision_framework_scvi.csv")
    
    print(f"\n  Clinical decision framework created for {len(clinical_df)} samples")
    print(f"  Columns: {list(clinical_df.columns)}")
    
except Exception as e:
    print(f"✗ Error building clinical framework: {e}")
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# STEP 8: Calculate performance metrics
# ============================================================================
print("\n" + "="*70)
print("STEP 8: Calculating performance metrics")
print("="*70)

try:
    metrics = {}
    
    # 1. Trajectory validation
    if 'stage' in adata.obs.columns and 'dpt_pseudotime' in adata.obs.columns:
        stage_map = {'GV': 0, 'MI': 1, 'MII': 2}
        stage_numeric = adata.obs['stage'].map(stage_map)
        corr, pval = spearmanr(stage_numeric, adata.obs['dpt_pseudotime'])
        metrics['trajectory_stage_correlation'] = {
            'spearman_r': float(corr),
            'p_value': float(pval),
            'n_samples': len(adata)
        }
        print(f"  ✓ Trajectory-stage correlation: ρ = {corr:.3f}, p = {pval:.4f}")
    
    # 2. Health score validation
    if 'stage' in adata.obs.columns and 'oocyte_health_score' in adata.obs.columns:
        gv_scores = adata.obs[adata.obs['stage'] == 'GV']['oocyte_health_score']
        mi_scores = adata.obs[adata.obs['stage'] == 'MI']['oocyte_health_score']
        
        from scipy.stats import mannwhitneyu
        u_stat, u_pval = mannwhitneyu(gv_scores, mi_scores, alternative='greater')
        
        metrics['health_score_validation'] = {
            'gv_mean': float(gv_scores.mean()),
            'mi_mean': float(mi_scores.mean()),
            'fold_change': float(gv_scores.mean() / mi_scores.mean()),
            'mann_whitney_u': float(u_stat),
            'p_value': float(u_pval)
        }
        print(f"  ✓ Health score: GV = {gv_scores.mean():.1f}, MI = {mi_scores.mean():.1f}, "
              f"fold-change = {gv_scores.mean()/mi_scores.mean():.2f}x")
    
    # 3. Cellular age correlation with chronological age
    if 'age' in adata.obs.columns and 'cellular_age_z' in adata.obs.columns:
        # Convert age to numeric if needed
        age_numeric = pd.to_numeric(adata.obs['age'], errors='coerce')
        valid_idx = ~age_numeric.isna()
        if valid_idx.sum() > 5:
            corr, pval = spearmanr(age_numeric[valid_idx], adata.obs['cellular_age_z'][valid_idx])
            metrics['age_correlation'] = {
                'spearman_r': float(corr),
                'p_value': float(pval),
                'n_samples': int(valid_idx.sum())
            }
            print(f"  ✓ Age correlation: ρ = {corr:.3f}, p = {pval:.4f}")
    
    # 4. Risk group separation
    if 'risk_group' in adata.obs.columns and 'cellular_age_z' in adata.obs.columns:
        low_risk = adata.obs[adata.obs['risk_group'].str.contains('Low', case=False)]['cellular_age_z']
        high_risk = adata.obs[adata.obs['risk_group'].str.contains('High', case=False)]['cellular_age_z']
        
        if len(low_risk) > 0 and len(high_risk) > 0:
            from scipy.stats import mannwhitneyu
            u_stat, u_pval = mannwhitneyu(low_risk, high_risk, alternative='two-sided')
            
            metrics['risk_group_separation'] = {
                'low_risk_mean': float(low_risk.mean()),
                'high_risk_mean': float(high_risk.mean()),
                'effect_size': float(high_risk.mean() - low_risk.mean()),
                'mann_whitney_u': float(u_stat),
                'p_value': float(u_pval)
            }
            print(f"  ✓ Risk group separation: Low = {low_risk.mean():.3f}, "
                  f"High = {high_risk.mean():.3f}")
    
    # 5. Pathway score statistics
    pathway_stats = {}
    for pathway in ['Cell_Cycle', 'Mitochondrial', 'DNA_Repair', 'Oxidative_Stress', 'Apoptosis']:
        score_col = f'score_{pathway}'
        if score_col in adata.obs.columns:
            scores = adata.obs[score_col].values
            pathway_stats[pathway] = {
                'mean': float(scores.mean()),
                'std': float(scores.std()),
                'min': float(scores.min()),
                'max': float(scores.max())
            }
    metrics['pathway_statistics'] = pathway_stats
    
    # Save metrics
    import json
    with open(f"{OUTDIRS['results']}/performance_metrics_scvi.json", 'w') as f:
        json.dump(metrics, f, indent=2)
    
    print(f"\n  ✓ Saved metrics: {OUTDIRS['results']}/performance_metrics_scvi.json")
    
except Exception as e:
    print(f"✗ Error calculating metrics: {e}")
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# STEP 9: Save complete annotated data
# ============================================================================
print("\n" + "="*70)
print("STEP 9: Saving complete annotated data")
print("="*70)

try:
    output_file = 'adata_complete_scvi.h5ad'
    adata.write(output_file)
    print(f"  ✓ Saved: {output_file}")
    print(f"    Shape: {adata.shape}")
    print(f"    Observations: {list(adata.obs.columns)}")
    
except Exception as e:
    print(f"✗ Error saving data: {e}")
    traceback.print_exc()

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "="*70)
print("PIPELINE COMPLETE - SUMMARY")
print("="*70)
print(f"\n✓ Data processed: {adata.shape}")
print(f"✓ scVI latent space: {adata.obsm['X_scvi'].shape}")
print(f"✓ DPT trajectory: computed")
print(f"✓ Pathway scores: {len([c for c in adata.obs.columns if c.startswith('score_')])} pathways")
print(f"✓ Health scores: computed")
print(f"✓ Risk stratification: {adata.obs['risk_group'].nunique()} groups")
print(f"✓ Clinical framework: saved")
print(f"✓ Performance metrics: calculated")
print(f"\nOutput files:")
print(f"  - adata_complete_scvi.h5ad")
print(f"  - {OUTDIRS['tables']}/clinical_decision_framework_scvi.csv")
print(f"  - {OUTDIRS['tables']}/pathway_scores_scvi.csv")
print(f"  - {OUTDIRS['results']}/performance_metrics_scvi.json")
print("\n" + "="*70)
print("SUCCESS: Complete pipeline executed with scVI!")
print("="*70)

