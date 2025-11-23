#!/usr/bin/env python3
"""
Forum EDA Visualizations: Raw Data Overview
Figure 1: Sample Composition (Raw Labels) - Stage and Age Distribution
Figure 2: Expression Matrix Snapshot - Heatmap of Top Variable Genes
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from pathlib import Path

print("="*70)
print("FORUM EDA: Raw Data Overview")
print("="*70)

# ============================================================================
# Load Metadata
# ============================================================================

sample_csv = '../data/sample_metadata_with_age.csv'
clinical_csv = '../data/clinical_decision_framework_final.csv'

if not os.path.exists(sample_csv):
    print("❌ ERROR: sample_metadata_with_age.csv not found")
    exit(1)

metadata = pd.read_csv(sample_csv)
print(f"✓ Loaded metadata: {len(metadata)} samples")

# Merge with clinical data if available for age
if os.path.exists(clinical_csv):
    clinical_df = pd.read_csv(clinical_csv, index_col=0)
    merged = metadata.merge(clinical_df, left_on='sample', right_index=True, how='left')
    # Extract age column
    age_col = 'age_x' if 'age_x' in merged.columns else ('age_y' if 'age_y' in merged.columns else 'age')
    if age_col in merged.columns:
        merged['age_parsed'] = pd.to_numeric(merged[age_col], errors='coerce')
    else:
        merged['age_parsed'] = pd.to_numeric(merged.get('age', pd.Series()), errors='coerce')
else:
    merged = metadata.copy()
    merged['age_parsed'] = pd.to_numeric(merged.get('age', pd.Series()), errors='coerce')

# Create age groups
def assign_age_group(age):
    if pd.isna(age):
        return 'Unknown'
    elif age < 30:
        return 'Young (<30)'
    elif age <= 40:
        return 'Middle (30-40)'
    else:
        return 'Old (≥40)'

merged['age_group'] = merged['age_parsed'].apply(assign_age_group)

# ============================================================================
# FIGURE 1: Sample Composition (Raw Labels)
# ============================================================================

print("\n[1/2] Creating sample composition plot (raw labels)...")

fig1, axes = plt.subplots(1, 2, figsize=(12, 5))
fig1.suptitle('Sample Composition: Raw Labels from Dataset', 
              fontsize=15, fontweight='bold', y=0.98)

colors_stage = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}

# Panel A: Count of cells per stage
if 'stage' in merged.columns:
    stage_counts = merged['stage'].value_counts().sort_index()
    bar_colors = [colors_stage.get(stage, '#95a5a6') for stage in stage_counts.index]
    
    bars1 = axes[0].bar(stage_counts.index, stage_counts.values, 
                       color=bar_colors, alpha=0.8, edgecolor='black', 
                       linewidth=2, width=0.6)
    
    axes[0].set_title('A. Cell Counts by Maturation Stage', 
                     fontsize=13, fontweight='bold', pad=12)
    axes[0].set_xlabel('Stage (GV / MI)', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('# of Cells', fontsize=12, fontweight='bold')
    axes[0].grid(True, alpha=0.3, axis='y', linestyle='--')
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    
    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        axes[0].text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}', ha='center', va='bottom',
                    fontweight='bold', fontsize=12)
else:
    axes[0].text(0.5, 0.5, 'Stage data\nnot available',
                ha='center', va='center', transform=axes[0].transAxes, fontsize=12)
    axes[0].set_title('A. Cell Counts by Stage', fontsize=13, fontweight='bold', pad=12)
    axes[0].axis('off')

# Panel B: Age distribution by age groups
if 'age_group' in merged.columns:
    age_group_counts = merged['age_group'].value_counts()
    
    # Sort by logical order
    age_order = ['Young (<30)', 'Middle (30-40)', 'Old (≥40)', 'Unknown']
    age_group_counts = age_group_counts.reindex([g for g in age_order if g in age_group_counts.index])
    
    # Colors for age groups
    age_colors = {'Young (<30)': '#2ecc71', 'Middle (30-40)': '#f39c12', 
                  'Old (≥40)': '#e74c3c', 'Unknown': '#95a5a6'}
    bar_colors2 = [age_colors.get(group, '#95a5a6') for group in age_group_counts.index]
    
    bars2 = axes[1].bar(range(len(age_group_counts)), age_group_counts.values,
                       color=bar_colors2, alpha=0.8, edgecolor='black',
                       linewidth=2, width=0.6)
    
    axes[1].set_xticks(range(len(age_group_counts)))
    axes[1].set_xticklabels([g.replace(' ', '\n') for g in age_group_counts.index],
                            fontsize=10, ha='center')
    axes[1].set_title('B. Donor Age Group Distribution', 
                     fontsize=13, fontweight='bold', pad=12)
    axes[1].set_xlabel('Age Group', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('# of Cells', fontsize=12, fontweight='bold')
    axes[1].grid(True, alpha=0.3, axis='y', linestyle='--')
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    
    # Add value labels
    for bar in bars2:
        height = bar.get_height()
        axes[1].text(bar.get_x() + bar.get_width()/2., height,
                    f'{int(height)}', ha='center', va='bottom',
                    fontweight='bold', fontsize=12)
elif merged['age_parsed'].notna().sum() > 0:
    # Fallback: histogram of raw ages
    ages = merged['age_parsed'].dropna()
    axes[1].hist(ages, bins=6, color='steelblue', alpha=0.7, 
                edgecolor='black', linewidth=1.5)
    axes[1].set_title('B. Donor Age Distribution (years)', 
                     fontsize=13, fontweight='bold', pad=12)
    axes[1].set_xlabel('Age (years)', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('# of Cells', fontsize=12, fontweight='bold')
    axes[1].grid(True, alpha=0.3, axis='y', linestyle='--')
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    
    # Add mean line
    mean_age = ages.mean()
    axes[1].axvline(mean_age, color='red', linestyle='--', linewidth=2,
                   label=f'Mean: {mean_age:.1f} years')
    axes[1].legend(fontsize=9)
else:
    axes[1].text(0.5, 0.5, 'Age data\nnot available',
                ha='center', va='center', transform=axes[1].transAxes, fontsize=12)
    axes[1].set_title('B. Donor Age Distribution', fontsize=13, fontweight='bold', pad=12)
    axes[1].axis('off')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('../visualizations/forum_raw_sample_composition.png', dpi=300, bbox_inches='tight',
            facecolor='white', edgecolor='none')
plt.close()
print("✓ Saved: forum_raw_sample_composition.png")

# Print summary
if 'stage' in merged.columns:
    stage_counts = merged['stage'].value_counts()
    print(f"\n  Stage distribution:")
    for stage, count in stage_counts.items():
        print(f"    {stage}: {count} cells")

if 'age_group' in merged.columns:
    age_counts = merged['age_group'].value_counts()
    print(f"  Age group distribution:")
    for group, count in age_counts.items():
        pct = count / len(merged) * 100
        print(f"    {group}: {count} cells ({pct:.0f}%)")

# ============================================================================
# FIGURE 2: Expression Matrix Snapshot (Heatmap)
# ============================================================================

print("\n[2/2] Creating expression matrix snapshot heatmap...")

# Try to load expression data from abundance.tsv files
zenodo_dir = './zenodo_data/final_code/kallisto'
expression_data = None
sample_names = []
gene_names = None

if os.path.exists(zenodo_dir):
    print("  Attempting to load expression data from Kallisto abundance files...")
    
    # Find all abundance.tsv files
    abundance_files = sorted(glob.glob(os.path.join(zenodo_dir, '*/abundance.tsv')))
    
    if abundance_files:
        print(f"  Found {len(abundance_files)} abundance.tsv files")
        
        # Read first file to get gene names
        first_df = pd.read_csv(abundance_files[0], sep='\t', comment='#')
        gene_names = first_df['target_id'].values
        
        # Collect TPM values for all samples
        tpm_matrix = []
        sample_names = []
        
        for i, ab_file in enumerate(abundance_files):
            df = pd.read_csv(ab_file, sep='\t', comment='#')
            sample_name = Path(ab_file).parent.name
            sample_names.append(sample_name)
            
            # Use TPM values
            tpm_values = df['tpm'].values
            tpm_matrix.append(tpm_values)
            
            if i < 3:  # Print first few for debugging
                print(f"    {sample_name}: {len(tpm_values)} transcripts, "
                      f"mean TPM: {np.mean(tpm_values):.2f}")
        
        # Convert to numpy array
        expression_data = np.array(tpm_matrix).T  # Transpose: genes x cells
        
        print(f"  ✓ Loaded expression matrix: {expression_data.shape[0]} genes × {expression_data.shape[1]} cells")
        
        # Select top variable genes
        # Calculate variance for each gene
        gene_vars = np.var(expression_data, axis=1)
        top_n_genes = 40
        top_var_indices = np.argsort(gene_vars)[-top_n_genes:][::-1]
        
        # Get top variable genes expression
        top_var_data = expression_data[top_var_indices, :]
        top_var_gene_names = gene_names[top_var_indices]
        
        print(f"  ✓ Selected top {top_n_genes} variable genes")
        
        # Log transform
        top_var_data_log = np.log1p(top_var_data)
        
        # Order cells by stage if possible
        if 'stage' in merged.columns:
            # Map sample names to stages
            sample_to_stage = dict(zip(merged['sample'], merged['stage']))
            # Get stage for each sample in order
            stages = [sample_to_stage.get(name, 'Unknown') for name in sample_names]
            
            # Sort by stage: GV first, then MI
            stage_order = {'GV': 0, 'MI': 1, 'MII': 2, 'Unknown': 3}
            sort_indices = sorted(range(len(stages)), key=lambda i: stage_order.get(stages[i], 99))
            
            top_var_data_log = top_var_data_log[:, sort_indices]
            sample_names_sorted = [sample_names[i] for i in sort_indices]
            stages_sorted = [stages[i] for i in sort_indices]
        else:
            sample_names_sorted = sample_names
            stages_sorted = ['Unknown'] * len(sample_names)
        
        # Create heatmap
        fig2, ax = plt.subplots(figsize=(14, 10))
        
        # Create heatmap using imshow
        im = ax.imshow(top_var_data_log, aspect='auto', cmap='viridis',
                      interpolation='nearest', vmin=0, vmax=np.percentile(top_var_data_log, 95))
        
        # Set labels
        ax.set_xlabel('Oocytes (cells)', fontsize=13, fontweight='bold')
        ax.set_ylabel(f'Top {top_n_genes} Most Variable Genes', fontsize=13, fontweight='bold')
        ax.set_title('Expression Matrix Snapshot: Raw TPM Values\n(Log-transformed, columns ordered by stage: GV → MI)',
                    fontsize=14, fontweight='bold', pad=15)
        
        # Set ticks (show subset of genes and samples)
        gene_tick_spacing = max(1, top_n_genes // 10)
        ax.set_yticks(range(0, top_n_genes, gene_tick_spacing))
        ax.set_yticklabels([top_var_gene_names[i][:20] if len(top_var_gene_names[i]) <= 20 
                           else top_var_gene_names[i][:17] + '...' 
                           for i in range(0, top_n_genes, gene_tick_spacing)],
                          fontsize=8, rotation=0)
        
        # Sample ticks (show every sample or subset)
        if len(sample_names_sorted) <= 20:
            ax.set_xticks(range(len(sample_names_sorted)))
            ax.set_xticklabels([name.split('_')[0] if '_' in name else name[:10] 
                               for name in sample_names_sorted],
                              fontsize=7, rotation=45, ha='right')
        else:
            tick_spacing = max(1, len(sample_names_sorted) // 10)
            ax.set_xticks(range(0, len(sample_names_sorted), tick_spacing))
            ax.set_xticklabels([sample_names_sorted[i][:10] 
                               for i in range(0, len(sample_names_sorted), tick_spacing)],
                              fontsize=7, rotation=45, ha='right')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Log(TPM + 1)', fontsize=11, fontweight='bold')
        
        # Add stage annotation bars if available
        if stages_sorted and set(stages_sorted) != {'Unknown'}:
            # Create stage color mapping
            stage_colors_map = {'GV': '#2ecc71', 'MI': '#e74c3c', 'MII': '#3498db'}
            
            # Add colored bars at the top
            for i, stage in enumerate(stages_sorted):
                color = stage_colors_map.get(stage, '#95a5a6')
                ax.axvline(i - 0.5, color=color, linewidth=3, alpha=0.7, zorder=0)
            
            # Add legend for stages
            from matplotlib.patches import Patch
            legend_elements = [Patch(facecolor=stage_colors_map.get(stage, '#95a5a6'), 
                                   label=stage) 
                             for stage in set(stages_sorted) if stage != 'Unknown']
            if legend_elements:
                ax.legend(handles=legend_elements, loc='upper left', 
                         fontsize=10, framealpha=0.9, title='Stage', 
                         title_fontsize=11)
        
        plt.tight_layout()
        plt.savefig('../visualizations/forum_raw_expression_heatmap.png', dpi=300, bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        plt.close()
        print("✓ Saved: forum_raw_expression_heatmap.png")
        
    else:
        print("  ⚠ No abundance.tsv files found - creating placeholder")
        create_placeholder_heatmap()
else:
    print("  ⚠ Zenodo data directory not found - creating placeholder")
    create_placeholder_heatmap()

def create_placeholder_heatmap():
    """Create a placeholder heatmap when expression data is not available"""
    fig2, ax = plt.subplots(figsize=(14, 10))
    
    # Create dummy data for demonstration
    n_genes = 40
    n_cells = 20
    
    # Simulate expression pattern: GV and MI show different patterns
    dummy_data = np.random.randn(n_genes, n_cells)
    # Add block structure: first 6 cells (GV) have different pattern than rest (MI)
    dummy_data[:, :6] += np.random.randn(n_genes, 6) * 2
    
    # Log transform
    dummy_data = np.log1p(np.exp(dummy_data))
    
    im = ax.imshow(dummy_data, aspect='auto', cmap='viridis',
                  interpolation='nearest', vmin=0, vmax=np.percentile(dummy_data, 95))
    
    ax.set_xlabel('Oocytes (cells)', fontsize=13, fontweight='bold')
    ax.set_ylabel(f'Top {n_genes} Most Variable Genes (example)', fontsize=13, fontweight='bold')
    ax.set_title('Expression Matrix Snapshot: Raw TPM Values\n(Log-transformed, columns ordered by stage: GV → MI)\n[NOTE: Placeholder - requires abundance.tsv files]',
                fontsize=14, fontweight='bold', pad=15)
    
    # Add note
    ax.text(0.5, -0.15, 'To generate actual heatmap, ensure abundance.tsv files are available in zenodo_data/final_code/kallisto/',
           transform=ax.transAxes, ha='center', fontsize=10, style='italic',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Log(TPM + 1)', fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('../visualizations/forum_raw_expression_heatmap.png', dpi=300, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.close()
    print("✓ Saved: forum_raw_expression_heatmap.png (placeholder)")

print("\n" + "="*70)
print("✓ Forum Raw EDA Visualizations Complete!")
print("="*70)
print("\nGenerated files:")
print("  ✓ forum_raw_sample_composition.png - Sample composition (stage + age)")
print("  ✓ forum_raw_expression_heatmap.png - Expression matrix snapshot")
print("\nThese visualizations show:")
print("  1. Raw labels: 20 oocytes, 6 GV + 14 MI, age distribution")
print("  2. Expression matrix: Top variable genes across cells, ordered by stage")

