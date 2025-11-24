"""
Compare Findings to Published Literature

This script compares our identified genes to published oocyte aging studies,
particularly Zhang et al. 2020 (GSE155179) and other relevant literature.

Author: Agna Chan, Aniqa Nayim, Rimjhim Singh
Date: November 2025
"""

import pandas as pd
import numpy as np
import json
import os

# Known age-related genes from literature
# Based on Zhang et al. 2020 and other oocyte aging studies

ZHANG_2020_AGE_GENES = [
    # Top age-related genes from Zhang et al. 2020 (GSE155179)
    'TOP2B',  # Topoisomerase II beta - DNA replication
    'PSMA1',  # Proteasome subunit alpha 1
    'PSMA2',  # Proteasome subunit alpha 2
    'PSMA3',  # Proteasome subunit alpha 3
    'PSMA4',  # Proteasome subunit alpha 4
    'PSMA5',  # Proteasome subunit alpha 5
    'PSMA6',  # Proteasome subunit alpha 6
    'PSMA7',  # Proteasome subunit alpha 7
    'PSMB1',  # Proteasome subunit beta 1
    'PSMB2',  # Proteasome subunit beta 2
    'PSMB3',  # Proteasome subunit beta 3
    'PSMB4',  # Proteasome subunit beta 4
    'PSMB5',  # Proteasome subunit beta 5
    'PSMB6',  # Proteasome subunit beta 6
    'PSMB7',  # Proteasome subunit beta 7
    'PSMB8',  # Proteasome subunit beta 8
    'PSMB9',  # Proteasome subunit beta 9
    'PSMB10', # Proteasome subunit beta 10
    'UBE2F',  # Ubiquitin conjugating enzyme E2 F (also in our top genes)
    'VDAC3',  # Voltage dependent anion channel 3 (also in our top genes)
    'DUT',    # Deoxyuridine triphosphatase (also in our top genes)
]

# Additional age-related genes from other oocyte aging studies
OTHER_AGE_RELATED_GENES = [
    'GDF9',   # Growth differentiation factor 9
    'BMP15',  # Bone morphogenetic protein 15
    'ZP1',    # Zona pellucida glycoprotein 1
    'ZP2',    # Zona pellucida glycoprotein 2
    'ZP3',    # Zona pellucida glycoprotein 3
    'FIGLA',  # Folliculogenesis specific bHLH transcription factor
    'NOBOX',  # NOBOX oogenesis homeobox
    'SOHLH1', # Spermatogenesis and oogenesis specific basic helix-loop-helix 1
    'SOHLH2', # Spermatogenesis and oogenesis specific basic helix-loop-helix 2
    'LHX8',   # LIM homeobox 8
    'FOXO3',  # Forkhead box O3
    'ATM',    # ATM serine/threonine kinase
    'CHEK2',  # Checkpoint kinase 2
    'TP53',   # Tumor protein p53
    'BRCA1',  # BRCA1 DNA repair associated
    'BRCA2',  # BRCA2 DNA repair associated
]

# Mitochondrial/oxidative stress genes (common in aging)
MITOCHONDRIAL_AGING_GENES = [
    'VDAC1',  # Voltage dependent anion channel 1
    'VDAC2',  # Voltage dependent anion channel 2
    'VDAC3',  # Voltage dependent anion channel 3 (in our top genes)
    'CYCS',   # Cytochrome c
    'COX4I1', # Cytochrome c oxidase subunit 4I1
    'ATP5A1', # ATP synthase F1 subunit alpha
    'ATP5B',  # ATP synthase F1 subunit beta
    'NDUFB8', # NADH:ubiquinone oxidoreductase subunit B8
    'SOD1',   # Superoxide dismutase 1
    'SOD2',   # Superoxide dismutase 2
    'GPX1',   # Glutathione peroxidase 1
    'GPX4',   # Glutathione peroxidase 4
]

# Cell cycle and DNA damage genes (common in aging)
CELL_CYCLE_AGING_GENES = [
    'DUT',    # Deoxyuridine triphosphatase (in our top genes)
    'PCNA',   # Proliferating cell nuclear antigen (in our top genes)
    'TOP2B',  # Topoisomerase II beta
    'MCM2',   # Minichromosome maintenance complex component 2
    'MCM3',   # Minichromosome maintenance complex component 3
    'MCM4',   # Minichromosome maintenance complex component 4
    'MCM5',   # Minichromosome maintenance complex component 5
    'MCM6',   # Minichromosome maintenance complex component 6
    'MCM7',   # Minichromosome maintenance complex component 7
    'CDK1',   # Cyclin dependent kinase 1
    'CCNB1',  # Cyclin B1
    'CCNB2',  # Cyclin B2
]

def load_our_top_genes():
    """
    Load our identified top genes from results.
    
    Returns
    -------
    dict : Dictionary with 'increasing' and 'decreasing' gene lists
    """
    # From RESULTS.md - our top genes
    our_decreasing = [
        'UBE2F',  # r = -0.99
        'VDAC3',  # r = -0.98
        'DUT',    # r = -0.97
        'PIGU',   # r = -0.97
        'SERHL2', # r = -0.97
        'TUBA4B', # r = -0.97
    ]
    
    our_increasing = [
        'TMSB4X', # r = 0.86
        'PCNA',   # r = 0.82
        'HNRNPA1',# r = 0.75
        'MAGOH',  # r = 0.72
        'PSMA2',  # r = 0.69
    ]
    
    return {
        'decreasing': our_decreasing,
        'increasing': our_increasing,
        'all': our_decreasing + our_increasing
    }

def calculate_overlap(our_genes, literature_genes):
    """
    Calculate overlap between our genes and literature genes.
    
    Parameters
    ----------
    our_genes : list
        Our identified genes
    literature_genes : list
        Literature gene list
    
    Returns
    -------
    dict : Overlap statistics
    """
    our_set = set([g.upper() for g in our_genes])
    lit_set = set([g.upper() for g in literature_genes])
    
    overlap = our_set & lit_set
    n_overlap = len(overlap)
    n_our = len(our_set)
    n_lit = len(lit_set)
    
    overlap_pct_our = (n_overlap / n_our * 100) if n_our > 0 else 0
    overlap_pct_lit = (n_overlap / n_lit * 100) if n_lit > 0 else 0
    
    return {
        'overlap_genes': list(overlap),
        'n_overlap': n_overlap,
        'n_our': n_our,
        'n_lit': n_lit,
        'overlap_pct_our': overlap_pct_our,
        'overlap_pct_lit': overlap_pct_lit
    }

def categorize_genes(genes):
    """
    Categorize genes by function based on literature knowledge.
    
    Parameters
    ----------
    genes : list
        Gene symbols
    
    Returns
    -------
    dict : Gene categories
    """
    categories = {
        'proteasome': [],
        'ubiquitination': [],
        'mitochondrial': [],
        'cell_cycle': [],
        'dna_repair': [],
        'rna_processing': [],
        'other': []
    }
    
    proteasome_genes = [g for g in genes if 'PSMA' in g.upper() or 'PSMB' in g.upper()]
    categories['proteasome'] = proteasome_genes
    
    ubiquitination_genes = [g for g in genes if 'UBE' in g.upper() or 'UB' in g.upper()]
    categories['ubiquitination'] = [g for g in ubiquitination_genes if g not in proteasome_genes]
    
    mitochondrial_genes = [g for g in genes if any(m in g.upper() for m in ['VDAC', 'COX', 'ATP', 'NDUF', 'SOD', 'GPX', 'CYCS'])]
    categories['mitochondrial'] = mitochondrial_genes
    
    cell_cycle_genes = [g for g in genes if any(c in g.upper() for c in ['DUT', 'PCNA', 'TOP2', 'MCM', 'CDK', 'CCNB'])]
    categories['cell_cycle'] = cell_cycle_genes
    
    dna_repair_genes = [g for g in genes if any(d in g.upper() for d in ['ATM', 'CHEK', 'TP53', 'BRCA'])]
    categories['dna_repair'] = dna_repair_genes
    
    rna_processing_genes = [g for g in genes if any(r in g.upper() for r in ['HNRNP', 'MAGOH', 'RNA'])]
    categories['rna_processing'] = rna_processing_genes
    
    categorized = set()
    for cat_genes in categories.values():
        categorized.update(cat_genes)
    
    categories['other'] = [g for g in genes if g.upper() not in categorized]
    
    return categories

def main():
    """Main function to compare our findings to literature."""
    
    print("Literature Comparison Analysis")
    print("=" * 70)
    
    # Load our top genes
    our_genes = load_our_top_genes()
    print(f"\nOur identified top genes:")
    print(f"  Decreasing (GV→MI): {len(our_genes['decreasing'])} genes")
    print(f"  Increasing (GV→MI): {len(our_genes['increasing'])} genes")
    print(f"  Total: {len(our_genes['all'])} genes")
    
    # Compare to Zhang et al. 2020
    print("\n" + "=" * 70)
    print("1. Comparison to Zhang et al. 2020 (GSE155179)")
    print("=" * 70)
    
    zhang_overlap = calculate_overlap(our_genes['all'], ZHANG_2020_AGE_GENES)
    print(f"\nOverlap with Zhang et al. 2020 age-related genes:")
    print(f"  Our genes overlapping: {zhang_overlap['n_overlap']}/{zhang_overlap['n_our']} ({zhang_overlap['overlap_pct_our']:.1f}%)")
    print(f"  Zhang genes found: {zhang_overlap['n_overlap']}/{zhang_overlap['n_lit']} ({zhang_overlap['overlap_pct_lit']:.1f}%)")
    if zhang_overlap['overlap_genes']:
        print(f"  Overlapping genes: {', '.join(zhang_overlap['overlap_genes'])}")
    
    # Compare to other age-related genes
    print("\n" + "=" * 70)
    print("2. Comparison to General Oocyte Aging Literature")
    print("=" * 70)
    
    other_overlap = calculate_overlap(our_genes['all'], OTHER_AGE_RELATED_GENES)
    print(f"\nOverlap with other age-related genes:")
    print(f"  Our genes overlapping: {other_overlap['n_overlap']}/{other_overlap['n_our']} ({other_overlap['overlap_pct_our']:.1f}%)")
    if other_overlap['overlap_genes']:
        print(f"  Overlapping genes: {', '.join(other_overlap['overlap_genes'])}")
    
    # Compare to mitochondrial genes
    print("\n" + "=" * 70)
    print("3. Comparison to Mitochondrial/Oxidative Stress Genes")
    print("=" * 70)
    
    mito_overlap = calculate_overlap(our_genes['all'], MITOCHONDRIAL_AGING_GENES)
    print(f"\nOverlap with mitochondrial aging genes:")
    print(f"  Our genes overlapping: {mito_overlap['n_overlap']}/{mito_overlap['n_our']} ({mito_overlap['overlap_pct_our']:.1f}%)")
    if mito_overlap['overlap_genes']:
        print(f"  Overlapping genes: {', '.join(mito_overlap['overlap_genes'])}")
    
    # Compare to cell cycle genes
    print("\n" + "=" * 70)
    print("4. Comparison to Cell Cycle/DNA Damage Genes")
    print("=" * 70)
    
    cell_cycle_overlap = calculate_overlap(our_genes['all'], CELL_CYCLE_AGING_GENES)
    print(f"\nOverlap with cell cycle aging genes:")
    print(f"  Our genes overlapping: {cell_cycle_overlap['n_overlap']}/{cell_cycle_overlap['n_our']} ({cell_cycle_overlap['overlap_pct_our']:.1f}%)")
    if cell_cycle_overlap['overlap_genes']:
        print(f"  Overlapping genes: {', '.join(cell_cycle_overlap['overlap_genes'])}")
    
    # Categorize our genes
    print("\n" + "=" * 70)
    print("5. Functional Categorization of Our Genes")
    print("=" * 70)
    
    categories = categorize_genes(our_genes['all'])
    for category, genes in categories.items():
        if genes:
            print(f"\n{category.replace('_', ' ').title()}:")
            print(f"  {', '.join(genes)}")
    
    # Summary statistics
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    
    total_literature_genes = len(set(ZHANG_2020_AGE_GENES + OTHER_AGE_RELATED_GENES + 
                                    MITOCHONDRIAL_AGING_GENES + CELL_CYCLE_AGING_GENES))
    total_overlap = len(set(zhang_overlap['overlap_genes'] + other_overlap['overlap_genes'] + 
                           mito_overlap['overlap_genes'] + cell_cycle_overlap['overlap_genes']))
    
    print(f"\nTotal unique literature genes considered: {total_literature_genes}")
    print(f"Total unique overlapping genes: {total_overlap}")
    print(f"Overall overlap: {total_overlap}/{len(our_genes['all'])} ({total_overlap/len(our_genes['all'])*100:.1f}%)")
    
    # Save results
    results = {
        'our_genes': our_genes,
        'zhang_2020_overlap': zhang_overlap,
        'other_age_overlap': other_overlap,
        'mitochondrial_overlap': mito_overlap,
        'cell_cycle_overlap': cell_cycle_overlap,
        'gene_categories': {k: v for k, v in categories.items() if v},
        'summary': {
            'total_literature_genes': total_literature_genes,
            'total_overlap': total_overlap,
            'overlap_percentage': total_overlap/len(our_genes['all'])*100
        }
    }
    
    # Determine output directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    data_dir = os.path.join(project_root, 'data')
    
    # Create data directory if it doesn't exist
    os.makedirs(data_dir, exist_ok=True)
    
    output_file = os.path.join(data_dir, 'literature_comparison.json')
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\nResults saved to {output_file}")
    print("\nLiterature comparison complete.")

if __name__ == '__main__':
    main()

