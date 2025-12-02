#!/usr/bin/env python3
"""
Add gene symbols to gene trajectory correlation results.
Maps transcript IDs (ENST...) to gene symbols using biomart file.
"""

import sys
import os
import pandas as pd
import numpy as np

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

print("Adding gene symbols to correlation results")

# Load biomart mapping file
biomart_path = os.path.join(BASE_DIR, 'zenodo_data/final_code/e106_biomart_HomoSapiens_Transcript_Gene_Names.txt')
print(f"\nLoading biomart file: {biomart_path}")

if not os.path.exists(biomart_path):
    print(f"Error: Biomart file not found at {biomart_path}")
    sys.exit(1)

biomart_df = pd.read_csv(biomart_path, sep='\t', header=0)
print(f"  Loaded {len(biomart_df)} transcript-to-gene mappings")

# Create mapping from transcript ID to gene symbol
# Transcript ID format: ENST00000555959.1
# Biomart has: "Transcript stable ID version" -> "Gene name"
transcript_to_gene = dict(zip(biomart_df['Transcript stable ID version'], biomart_df['Gene name']))

# Also create mapping without version number (e.g., ENST00000555959 -> Gene name)
# Some transcript IDs might not have version suffixes
transcript_base_to_gene = {}
for transcript_id, gene_name in transcript_to_gene.items():
    # Remove version suffix
    base_id = transcript_id.split('.')[0]
    if base_id not in transcript_base_to_gene:
        transcript_base_to_gene[base_id] = gene_name

print(f"  Created mapping for {len(transcript_to_gene)} transcripts")
print(f"  Created base ID mapping for {len(transcript_base_to_gene)} transcripts")

# Load correlation results
corr_path = os.path.join(BASE_DIR, 'pipeline_results_scvi/tables/gene_trajectory_correlations.csv')
print(f"\nLoading correlation results: {corr_path}")

if not os.path.exists(corr_path):
    print(f"Error: Correlation file not found at {corr_path}")
    sys.exit(1)

corr_df = pd.read_csv(corr_path)
print(f"  Loaded {len(corr_df)} gene correlations")

# Map transcript IDs to gene symbols
print("\nMapping transcript IDs to gene symbols...")

def get_gene_symbol(transcript_id):
    """Get gene symbol for a transcript ID."""
    # Try exact match first
    if transcript_id in transcript_to_gene:
        return transcript_to_gene[transcript_id]
    
    # Try without version suffix
    base_id = transcript_id.split('.')[0]
    if base_id in transcript_base_to_gene:
        return transcript_base_to_gene[base_id]
    
    # If no match, return transcript ID itself
    return transcript_id

corr_df['gene_symbol'] = corr_df['gene'].apply(get_gene_symbol)

# Count how many were successfully mapped
mapped_count = (corr_df['gene_symbol'] != corr_df['gene']).sum()
print(f"  Successfully mapped {mapped_count}/{len(corr_df)} transcripts to gene symbols")

# Reorder columns: gene_symbol first, then transcript_id (gene), then others
cols = ['gene_symbol', 'gene'] + [c for c in corr_df.columns if c not in ['gene_symbol', 'gene']]
corr_df = corr_df[cols]

# Rename 'gene' to 'transcript_id' for clarity
corr_df = corr_df.rename(columns={'gene': 'transcript_id'})

# Sort by absolute correlation (highest first)
corr_df = corr_df.sort_values('abs_correlation', ascending=False).reset_index(drop=True)

# Save updated results
output_path = os.path.join(BASE_DIR, 'pipeline_results_scvi/tables/gene_trajectory_correlations_with_symbols.csv')
corr_df.to_csv(output_path, index=False)
print(f"Saved updated correlation results to: {output_path}")

# Show top genes
print("\nTop 10 correlated genes (with gene symbols)")

print("\nTop Decreasing Genes (GV→MI):")
top_decreasing = corr_df[corr_df['correlation'] < 0].head(10)
for idx, row in top_decreasing.iterrows():
    print(f"  {row['gene_symbol']:15s} ({row['transcript_id']:20s}): ρ = {row['correlation']:7.3f}, FDR = {row['fdr']:.4f}")

print("\nTop Increasing Genes (GV→MI):")
top_increasing = corr_df[corr_df['correlation'] > 0].head(10)
for idx, row in top_increasing.iterrows():
    print(f"  {row['gene_symbol']:15s} ({row['transcript_id']:20s}): ρ = {row['correlation']:7.3f}, FDR = {row['fdr']:.4f}")

# Summary statistics by gene symbol (in case multiple transcripts map to same gene)
print("\n" + "="*70)
print("SUMMARY BY GENE SYMBOL")
print("="*70)

# Group by gene symbol and take the transcript with highest absolute correlation
def get_max_corr(x):
    """Get correlation value for transcript with maximum absolute correlation."""
    idx = x.abs().idxmax()
    return x.loc[idx]

def get_max_corr_transcript(x):
    """Get transcript ID for transcript with maximum absolute correlation."""
    idx = x.abs().idxmax()
    return corr_df.loc[idx, 'transcript_id']

gene_summary = corr_df.groupby('gene_symbol').agg({
    'correlation': get_max_corr,  # Correlation of most correlated transcript
    'abs_correlation': 'max',  # Maximum absolute correlation
    'fdr': 'min',  # Minimum FDR (most significant)
}).reset_index()

# Get transcript ID for the most correlated transcript per gene
gene_summary['transcript_id'] = corr_df.groupby('gene_symbol')['abs_correlation'].apply(
    lambda x: corr_df.loc[x.idxmax(), 'transcript_id']
).reset_index(drop=True)

gene_summary = gene_summary.sort_values('abs_correlation', ascending=False)

# Count unique genes with significant correlations
n_significant_genes = ((gene_summary['abs_correlation'] > 0.7) & (gene_summary['fdr'] < 0.1)).sum()
print(f"\nUnique genes with |ρ| > 0.7 and FDR < 0.1: {n_significant_genes}")

print(f"\nTop 20 genes by absolute correlation:")
print(gene_summary.head(20)[['gene_symbol', 'correlation', 'fdr']].to_string(index=False))

# Save gene-level summary
gene_summary_path = os.path.join(BASE_DIR, 'pipeline_results_scvi/tables/gene_symbol_correlations_summary.csv')
gene_summary.to_csv(gene_summary_path, index=False)
print(f"Saved gene-level summary to: {gene_summary_path}")

print("\nAnalysis complete")

