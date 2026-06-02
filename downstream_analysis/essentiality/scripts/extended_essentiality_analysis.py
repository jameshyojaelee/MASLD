#!/usr/bin/env python3
"""
Extended Essentiality Analysis
==============================

This script generates publication-ready essentiality visualizations for Core DEGs.
It uses both Liang Cas13 screen data and DepMap Cas9 data.

Generated Plots:
1. Essentiality vs LFC Volcano Plot
2. Human-Mouse Concordance for Essential Genes
3. Lineage-Specific Essentiality Heatmap (Liver cell lines)

Author: Sanjana Lab Computational Pipeline
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
from pathlib import Path
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

# ============== CONFIGURATION ==============
ESSENTIALITY_DIR = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/essentiality")
OUTPUT_DIR = ESSENTIALITY_DIR / "plots"


# Files
CORE_DEGS_FILE = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/final_core_degs.csv")
DEPMAP_FILE = ESSENTIALITY_DIR / "CRISPRGeneEffect.csv"
MODEL_FILE = ESSENTIALITY_DIR / "Model.csv"
UNIFIED_DATA = OUTPUT_DIR / "Unified_Essentiality_Combined.csv"

# Thresholds
ESSENTIAL_THRESH = -0.5  # DepMap convention

# Cell lines to exclude from analysis
EXCLUDED_CELL_LINES = ['KMCH1', 'SKHEP1']

# Colors - Sanjana Lab palette
COLORS = {
    'primary': '#EC407A',      # Pink
    'secondary': '#FF7043',    # Orange
    'tertiary': '#7E57C2',     # Purple
    'background': '#9E9E9E',   # Gray
    'essential': '#D81B60',    # Deep pink
    'non_essential': '#42A5F5' # Blue
}

sns.set_context("paper", font_scale=1.8)
sns.set_style("whitegrid", {'axes.grid': False})

# Configure fonts for Illustrator (Type 42 = TrueType)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Font Family (Helvetica/Arial)
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'sans-serif']

# Strict Font Sizes (Journal Standards)
matplotlib.rcParams['font.size'] = 7           # Default text size
matplotlib.rcParams['axes.titlesize'] = 8      # Title size
matplotlib.rcParams['axes.labelsize'] = 8      # X and Y label size
matplotlib.rcParams['xtick.labelsize'] = 6     # X tick size
matplotlib.rcParams['ytick.labelsize'] = 6     # Y tick size
matplotlib.rcParams['legend.fontsize'] = 6     # Legend size

# ============== DATA LOADING ==============
print("Loading data...")

# Load unified essentiality data (previously computed)
if UNIFIED_DATA.exists():
    unified_df = pd.read_csv(UNIFIED_DATA)
    print(f"  Loaded unified data: {len(unified_df)} genes")
else:
    print("  Unified data not found, loading from source...")
    unified_df = pd.DataFrame()

# Load Core DEGs
core_df = pd.read_csv(CORE_DEGS_FILE)
print(f"  Loaded Core DEGs: {len(core_df)} entries")

# Load DepMap data
print("  Loading DepMap CRISPRGeneEffect (this may take a moment)...")
depmap_df = pd.read_csv(DEPMAP_FILE, index_col=0)
print(f"  Loaded DepMap: {depmap_df.shape[0]} cell lines × {depmap_df.shape[1]} genes")

# Load Model metadata
model_df = pd.read_csv(MODEL_FILE)
print(f"  Loaded Model metadata: {len(model_df)} models")

# ============== PREPROCESSING ==============
# Extract gene symbols from column names (format: "GENE (12345)")
def extract_symbol(col):
    if ' (' in col:
        return col.split(' (')[0]
    return col

depmap_df.columns = [extract_symbol(c) for c in depmap_df.columns]

# Mean essentiality per gene across all cell lines
mean_essentiality = depmap_df.mean(axis=0).reset_index()
mean_essentiality.columns = ['Gene', 'Mean_Essentiality']
mean_essentiality = mean_essentiality.sort_values('Mean_Essentiality')
mean_essentiality['Rank'] = range(1, len(mean_essentiality) + 1)
print(f"  Computed mean essentiality for {len(mean_essentiality)} genes")

# Core gene set (human + mapped mouse)
human_genes = core_df[core_df['species'] == 'human']['gene_symbol'].dropna().unique()
# For mouse, we need orthologs - use the 'human_ortholog_gene' if available
if 'human_ortholog_gene' in core_df.columns:
    mouse_orthologs = core_df[core_df['species'] == 'mouse']['human_ortholog_gene'].dropna().unique()
    core_symbols = set(human_genes) | set(mouse_orthologs)
else:
    core_symbols = set(human_genes)
print(f"  Core DEG gene symbols: {len(core_symbols)}")

# Add Core DEG flag
mean_essentiality['IsCoreDEG'] = mean_essentiality['Gene'].isin(core_symbols)
core_in_depmap = mean_essentiality[mean_essentiality['IsCoreDEG']]
print(f"  Core DEGs found in DepMap: {len(core_in_depmap)}")

# ============== PLOT 1: ESSENTIALITY VS LFC VOLCANO ==============
print("\nPlot 1: Essentiality vs LFC Volcano...")

# Need to merge with LFC data from Core DEGs
if 'log2fc_mean' in core_df.columns or 'log2FoldChange' in core_df.columns:
    lfc_col = 'log2fc_mean' if 'log2fc_mean' in core_df.columns else 'log2FoldChange'
    
    # Get human Core DEGs with LFC
    human_core = core_df[core_df['species'] == 'human'][['gene_symbol', lfc_col]].dropna()
    human_core = human_core.drop_duplicates(subset='gene_symbol')
    
    # Merge with essentiality
    volcano_df = human_core.merge(mean_essentiality[['Gene', 'Mean_Essentiality']], 
                                   left_on='gene_symbol', right_on='Gene', how='inner')
    
    if len(volcano_df) > 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Categorize
        volcano_df['Category'] = 'Other'
        volcano_df.loc[(volcano_df['Mean_Essentiality'] < ESSENTIAL_THRESH) & 
                       (volcano_df[lfc_col] > 0), 'Category'] = 'Essential + Upregulated'
        volcano_df.loc[(volcano_df['Mean_Essentiality'] < ESSENTIAL_THRESH) & 
                       (volcano_df[lfc_col] < 0), 'Category'] = 'Essential + Downregulated'
        volcano_df.loc[(volcano_df['Mean_Essentiality'] >= ESSENTIAL_THRESH) & 
                       (abs(volcano_df[lfc_col]) > 1), 'Category'] = 'Non-Essential + High LFC'
        
        color_map = {
            'Essential + Upregulated': COLORS['primary'],
            'Essential + Downregulated': COLORS['tertiary'],
            'Non-Essential + High LFC': COLORS['secondary'],
            'Other': COLORS['background']
        }
        
        for cat, color in color_map.items():
            subset = volcano_df[volcano_df['Category'] == cat]
            ax.scatter(subset[lfc_col], subset['Mean_Essentiality'], 
                      c=color, label=f'{cat} (n={len(subset)})', alpha=0.7, s=40)
        
        ax.axhline(ESSENTIAL_THRESH, color='gray', linestyle='--', alpha=0.5)
        ax.axvline(0, color='gray', linestyle='-', alpha=0.3)
        
        ax.set_xlabel('Log2 Fold Change (Disease vs Control)', fontsize=14)
        ax.set_ylabel('Essentiality Score (DepMap)', fontsize=14)
        ax.set_title('Essentiality vs Differential Expression', fontsize=16, fontweight='bold')
        ax.legend(loc='lower left', fontsize=10, frameon=False)
        sns.despine()
        plt.tight_layout()
        # PNG removed
        plt.savefig(OUTPUT_DIR / "New_Essentiality_vs_LFC_Volcano.pdf", bbox_inches='tight')
        plt.close()
        print("  Saved: New_Essentiality_vs_LFC_Volcano.pdf")
    else:
        print("  Skipped: No overlap between Core DEGs and DepMap")
else:
    print("  Skipped: No LFC column in Core DEGs")

# ============== PLOT 2: HUMAN-MOUSE CONCORDANCE ==============
print("\nPlot 2: Human-Mouse Essentiality Concordance...")

# Check if we have both human and mouse data with orthologs
if 'species' in core_df.columns:
    human_core_genes = core_df[core_df['species'] == 'human']['gene_symbol'].dropna().unique()
    mouse_core_genes = core_df[core_df['species'] == 'mouse']
    
    if 'human_ortholog_gene' in core_df.columns:
        mouse_orthologs = core_df[core_df['species'] == 'mouse'][['gene_symbol', 'human_ortholog_gene']].dropna()
        
        # Get essentiality for human genes and their mouse orthologs
        human_ess = mean_essentiality[mean_essentiality['Gene'].isin(human_core_genes)].copy()
        human_ess = human_ess.rename(columns={'Gene': 'Human_Gene', 'Mean_Essentiality': 'Human_Essentiality'})
        
        ortho_ess = mean_essentiality[mean_essentiality['Gene'].isin(mouse_orthologs['human_ortholog_gene'])].copy()
        ortho_ess = ortho_ess.rename(columns={'Gene': 'Human_Gene', 'Mean_Essentiality': 'Mouse_Ortholog_Essentiality'})
        
        concordance_df = human_ess[['Human_Gene', 'Human_Essentiality']].merge(
            ortho_ess[['Human_Gene', 'Mouse_Ortholog_Essentiality']], on='Human_Gene', how='inner')
        
        if len(concordance_df) > 10:
            fig, ax = plt.subplots(figsize=(9, 9))
            
            ax.scatter(concordance_df['Human_Essentiality'], 
                      concordance_df['Mouse_Ortholog_Essentiality'],
                      c=COLORS['primary'], alpha=0.6, s=50)
            
            # Add regression line
            from scipy import stats
            slope, intercept, r, p, se = stats.linregress(concordance_df['Human_Essentiality'],
                                                           concordance_df['Mouse_Ortholog_Essentiality'])
            x_line = np.linspace(concordance_df['Human_Essentiality'].min(), 
                                concordance_df['Human_Essentiality'].max(), 100)
            ax.plot(x_line, slope * x_line + intercept, 
                   color=COLORS['secondary'], linewidth=2, linestyle='--')
            
            ax.axhline(ESSENTIAL_THRESH, color='gray', linestyle=':', alpha=0.5)
            ax.axvline(ESSENTIAL_THRESH, color='gray', linestyle=':', alpha=0.5)
            
            ax.set_xlabel('Human DEG Essentiality', fontsize=14)
            ax.set_ylabel('Mouse Ortholog Essentiality', fontsize=14)
            ax.set_title('Human-Mouse Essentiality Concordance', fontsize=16, fontweight='bold')
            ax.text(0.05, 0.95, f'r = {r:.3f}\np = {p:.2e}\nn = {len(concordance_df)}',
                   transform=ax.transAxes, fontsize=11, va='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            sns.despine()
            plt.tight_layout()
            # PNG removed
            plt.savefig(OUTPUT_DIR / "New_Human_Mouse_Concordance.pdf", bbox_inches='tight')
            plt.close()
            print(f"  Saved: New_Human_Mouse_Concordance.pdf (r={r:.3f})")
        else:
            print("  Skipped: Insufficient concordance data")
    else:
        print("  Skipped: No human_ortholog_gene column")
else:
    print("  Skipped: No species column")

# ============== PLOT 3: LIVER LINEAGE HEATMAP ==============
print("\nPlot 3: Liver-Specific Essentiality Heatmap...")

# File path for Quadrants data
QUADRANTS_FILE = OUTPUT_DIR / "Essentiality_Quadrants_All.csv"

# Filter to liver cell lines
liver_models = model_df[model_df['OncotreeLineage'].str.contains('Liver', case=False, na=False)]

# EXCLUDE specified cell lines
excluded_mask = liver_models['StrippedCellLineName'].isin(EXCLUDED_CELL_LINES)
n_excluded = excluded_mask.sum()
liver_models = liver_models[~excluded_mask]
print(f"  Excluded {n_excluded} cell lines: {EXCLUDED_CELL_LINES}")

liver_model_ids = liver_models['ModelID'].values

# Get liver cell lines present in DepMap
liver_in_depmap = [m for m in liver_model_ids if m in depmap_df.index]
print(f"  Found {len(liver_in_depmap)} liver cell lines in DepMap (after exclusion)")

if len(liver_in_depmap) > 0:
    # Load target genes from Quadrants file if it exists, otherwise fall back to Core DEGs
    target_genes = []
    if QUADRANTS_FILE.exists():
        print(f"  Loading target genes from: {QUADRANTS_FILE}")
        quad_df = pd.read_csv(QUADRANTS_FILE)
        if 'Gene' in quad_df.columns:
            target_genes = quad_df['Gene'].unique()
            print(f"  Found {len(target_genes)} genes in Quadrants file")
        else:
            print("  Warning: 'Gene' column not found in Quadrants file. Using essential Core DEGs.")
    else:
        print(f"  Quadrants file not found at {QUADRANTS_FILE}")

    # Fallback if no genes loaded
    if len(target_genes) == 0:
        target_genes = core_in_depmap[core_in_depmap['Mean_Essentiality'] < ESSENTIAL_THRESH]['Gene'].values

    # Filter for genes present in DepMap
    valid_genes = [g for g in target_genes if g in depmap_df.columns]
    
    # Sort genes by mean essentiality in liver lines for better visualization
    liver_data_raw = depmap_df.loc[liver_in_depmap, valid_genes]
    gene_means = liver_data_raw.mean(axis=0).sort_values()
    sorted_genes = gene_means.index.tolist()
    
    print(f"  Plotting {len(sorted_genes)} genes (validated in DepMap)")
    
    if len(sorted_genes) > 5:
        heatmap_data = depmap_df.loc[liver_in_depmap, sorted_genes]
        
        # Add cell line names
        liver_names = liver_models.set_index('ModelID').loc[liver_in_depmap, 'StrippedCellLineName'].values
        heatmap_data.index = liver_names
        
        # Dynamic Figure Size
        fig_width = max(12, len(sorted_genes) * 0.15)
        fig_height = max(8, len(liver_in_depmap) * 0.4)
        
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        # Create heatmap
        cmap = sns.color_palette("RdPu_r", as_cmap=True)
        sns.heatmap(heatmap_data, cmap=cmap, center=0, 
                   xticklabels=True, yticklabels=True,
                   cbar_kws={'label': 'Essentiality Score', 'shrink': 0.5})
        
        ax.set_xlabel('Genes (Sorted by Mean Liver Essentiality)', fontsize=14)
        ax.set_ylabel('Liver Cell Lines', fontsize=14)
        ax.set_title(f'Liver-Specific Essentiality ({len(sorted_genes)} Genes)', fontsize=16, fontweight='bold')
        
        # Adjust tick font size based on number of genes
        tick_fontsize = 4 if len(sorted_genes) > 200 else 9
        plt.xticks(rotation=90, ha='center', fontsize=tick_fontsize)
        plt.yticks(fontsize=10)
        
        plt.tight_layout()
        
        # Save
        # PNG removed
        plt.savefig(OUTPUT_DIR / "New_Liver_Essentiality_Heatmap.pdf", bbox_inches='tight')
        plt.close()
        print(f"  Saved: New_Liver_Essentiality_Heatmap.pdf ({len(liver_in_depmap)} lines × {len(sorted_genes)} genes)")
    else:
        print("  Skipped: Insufficient genes for heatmap")
else:
    print("  Skipped: Insufficient liver cell lines")

# ============== SUMMARY ==============
print("\n" + "="*60)
print("Extended Essentiality Analysis Complete!")
print("="*60)
print(f"\nPlots saved to: {OUTPUT_DIR}")
print("\nGenerated:")
print("  1. New_Essentiality_vs_LFC_Volcano.pdf")
print("  2. New_Human_Mouse_Concordance.pdf")
print("  3. New_Liver_Essentiality_Heatmap.pdf")
print(f"\nExcluded cell lines: {EXCLUDED_CELL_LINES}")
