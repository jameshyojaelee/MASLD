#!/usr/bin/env python3
"""
Extended Deconvolution Analysis: 5 New Creative Visualizations
================================================================

This script generates 5 new publication-ready plots that complement the existing
MuSiC deconvolution analyses. Each visualization addresses a unique biological question.

New Plots:
1. Disease Progression Trajectory (cell type changes across fibrosis stages)
2. Human vs Mouse Cell Type Comparison (species concordance)
3. Cell Type Correlation Network (which cell types co-vary)
4. Individual Sample Variability (raincloud plots per cell type)
5. Top Variable Cell Types Radar Chart

Author: Sanjana Lab Computational Pipeline
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# ============== CONFIGURATION ==============
BASE_DIR = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/deconvolution")
RESULTS_DIR = BASE_DIR / "results"
OUTPUT_DIR = RESULTS_DIR / "summary_plots"
METADATA_DIR = BASE_DIR / "metadata"

# Datasets
DATASETS = {
    "GSE130970": "Human (Hoang)",
    "GSE135251": "Human (Govaere)",
    "GSE156918": "Mouse MCD (Paquette)",
    "GSE205974": "Mouse MCD (Yue)",
    "inhouse_MCD": "Cas13 mouse MCD"
}

import argparse

parser = argparse.ArgumentParser(description='Generate Extended Deconvolution Plots')
parser.add_argument('--method', type=str, default='music', choices=['music', 'bayesprism'],
                    help='Deconvolution method to analyze')
args = parser.parse_args()
METHOD = args.method

# Colors - Sanjana Lab palette
COLORS = {
    'primary': '#EC407A',      # Pink
    'secondary': '#FF7043',    # Orange
    'tertiary': '#7E57C2',     # Purple
    'quaternary': '#26A69A',   # Teal
    'background': '#9E9E9E',   # Gray
}

# Cell type color palette
CELL_COLORS = sns.color_palette("husl", 15)

# Set context and style
sns.set_context("paper") 
sns.set_style("whitegrid", {'axes.grid': False})

# CRITICAL: Configure fonts for Illustrator (Type 42 = TrueType)
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Font Family
# Note: 'Helvetica' removed from runtime list to prevent Linux segfaults
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']

# Strict Font Sizes (Journal Standards)
plt.rcParams['font.size'] = 7           # Default text size
plt.rcParams['axes.titlesize'] = 8      # Title size
plt.rcParams['axes.labelsize'] = 8      # X and Y label size
plt.rcParams['xtick.labelsize'] = 6     # X tick size
plt.rcParams['ytick.labelsize'] = 6     # Y tick size
plt.rcParams['legend.fontsize'] = 6     # Legend size

# ============== DATA LOADING ==============
print("Loading deconvolution results...")

all_data = []
for dataset in DATASETS.keys():
    # Load based on method
    if METHOD == 'music':
        weighted_file = RESULTS_DIR / dataset / f"{dataset}_music_prop_weighted.tsv"
        all_file = RESULTS_DIR / dataset / f"{dataset}_music_prop_all.tsv"
        
        if weighted_file.exists():
            df = pd.read_csv(weighted_file, sep="\t", index_col=0)
        elif all_file.exists():
            df = pd.read_csv(all_file, sep="\t", index_col=0)
        else:
            print(f"  Warning: No results found for {dataset}")
            continue
    elif METHOD == 'bayesprism':
        file_path = RESULTS_DIR / dataset / f"{dataset}_bayesprism_proportions.tsv"
        if file_path.exists():
            df = pd.read_csv(file_path, sep="\t", index_col=0)
        else:
            print(f"  Warning: No results found for {dataset}")
            continue
    
    # Melt to long format
    df_long = df.reset_index().melt(id_vars=[df.index.name or 'index'], 
                                     var_name='CellType', value_name='Fraction')
    df_long = df_long.rename(columns={df.index.name or 'index': 'Sample'})
    df_long['Dataset'] = dataset
    df_long['Species'] = 'Human' if 'GSE130970' in dataset or 'GSE135251' in dataset else 'Mouse'
    all_data.append(df_long)
    print(f"  Loaded {dataset}: {len(df)} samples × {df.shape[1]} cell types")

if not all_data:
    print("ERROR: No data loaded!")
    exit(1)

combined_df = pd.concat(all_data, ignore_index=True)
print(f"\nCombined data: {len(combined_df)} rows")

# Get unique cell types
cell_types = combined_df['CellType'].unique()
print(f"Cell types: {len(cell_types)}")

# ============== LOAD METADATA FOR DISEASE STAGING ==============
# Try to load metadata for fibrosis staging
metadata_files = {
    "GSE130970": BASE_DIR.parent / "patient_RNAseq" / "data" / "metadata" / "GSE130970_SraRunTable.csv",
    "GSE135251": BASE_DIR.parent / "patient_RNAseq" / "data" / "metadata" / "GSE135251_SraRunTable.csv",
}

staging_data = {}
for dataset, meta_path in metadata_files.items():
    if meta_path.exists():
        meta = pd.read_csv(meta_path)
        # Find fibrosis column
        fib_col = [c for c in meta.columns if 'fibrosis' in c.lower()]
        sample_col = [c for c in meta.columns if c.lower() in ['run', 'sample', 'sample_name']]
        if fib_col and sample_col:
            staging_data[dataset] = meta[[sample_col[0], fib_col[0]]].rename(
                columns={sample_col[0]: 'Sample', fib_col[0]: 'Fibrosis'})
            print(f"  Loaded staging for {dataset}")

# ============== PLOT 1: DISEASE PROGRESSION TRAJECTORY ==============
print("\nPlot 1: Disease Progression Trajectory...")

if staging_data:
    # Find datasets with staging
    staged_datasets = list(staging_data.keys())
    if staged_datasets:
        fig, axes = plt.subplots(1, len(staged_datasets), figsize=(5*len(staged_datasets), 4.5))
        if len(staged_datasets) == 1:
            axes = [axes]
        
        for ax, dataset in zip(axes, staged_datasets):
            # Get deconvolution data
            data = combined_df[combined_df['Dataset'] == dataset].copy()
            staging = staging_data[dataset]
            
            # Merge with staging
            data = data.merge(staging, on='Sample', how='left')
            data = data.dropna(subset=['Fibrosis'])
            
            if len(data) > 0:
                # Convert fibrosis to numeric if possible
                data['Fibrosis'] = pd.to_numeric(data['Fibrosis'], errors='coerce')
                data = data.dropna(subset=['Fibrosis'])
                
                # Get top 5 cell types by variance
                ct_var = data.groupby('CellType')['Fraction'].std().nlargest(5).index
                data_top = data[data['CellType'].isin(ct_var)]
                
                # Plot trajectory
                for i, ct in enumerate(ct_var):
                    ct_data = data_top[data_top['CellType'] == ct]
                    mean_by_fib = ct_data.groupby('Fibrosis')['Fraction'].mean()
                    sem_by_fib = ct_data.groupby('Fibrosis')['Fraction'].sem()
                    
                    ax.errorbar(mean_by_fib.index, mean_by_fib.values, 
                               yerr=sem_by_fib.values, 
                               label=ct, marker='o', markersize=8,
                               linewidth=2, capsize=3, color=CELL_COLORS[i])
                
                ax.set_xlabel('Fibrosis Stage', fontsize=12)
                ax.set_ylabel('Cell Type Fraction', fontsize=12)
                ax.set_title(f'{DATASETS[dataset]}', fontsize=14, fontweight='bold')
                ax.legend(loc='best', fontsize=9, frameon=False)
        
        plt.suptitle(f'Cell Type Changes Across Fibrosis Progression ({METHOD})', 
                    fontsize=16, fontweight='bold', y=1.02)
        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / f"{METHOD}_New_Disease_Progression_Trajectory.pdf", bbox_inches='tight')
        plt.close()
        print("  Saved: New_Disease_Progression_Trajectory.pdf")
else:
    print("  Skipped: No staging metadata available")

# ============== PLOT 2: HUMAN VS MOUSE COMPARISON ==============
print("\nPlot 2: Human vs Mouse Cell Type Comparison...")

# Calculate mean fractions per species per cell type
species_means = combined_df.groupby(['Species', 'CellType'])['Fraction'].mean().reset_index()
species_pivot = species_means.pivot(index='CellType', columns='Species', values='Fraction').fillna(0)

if 'Human' in species_pivot.columns and 'Mouse' in species_pivot.columns:
    fig, ax = plt.subplots(figsize=(6, 6))
    
    ax.scatter(species_pivot['Human'], species_pivot['Mouse'], 
              s=100, c=COLORS['primary'], alpha=0.7, edgecolors='white', linewidth=1)
    
    # Add labels for top cell types
    for ct in species_pivot.index:
        if species_pivot.loc[ct, 'Human'] > 0.05 or species_pivot.loc[ct, 'Mouse'] > 0.05:
            ax.annotate(ct, (species_pivot.loc[ct, 'Human'], species_pivot.loc[ct, 'Mouse']),
                       fontsize=9, ha='left', va='bottom')
    
    # Add diagonal line
    max_val = max(species_pivot['Human'].max(), species_pivot['Mouse'].max())
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, linewidth=1)
    
    # Correlation
    from scipy import stats
    r, p = stats.pearsonr(species_pivot['Human'], species_pivot['Mouse'])
    ax.text(0.05, 0.95, f'r = {r:.3f}\np = {p:.2e}',
           transform=ax.transAxes, fontsize=12, va='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_xlabel('Human NAFLD (mean fraction)', fontsize=14)
    ax.set_ylabel('Mouse MCD (mean fraction)', fontsize=14)
    ax.set_title(f'Human vs Mouse Cell Type Concordance ({METHOD})', fontsize=16, fontweight='bold')
    sns.despine()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / f"{METHOD}_New_Human_Mouse_Comparison.pdf", bbox_inches='tight')
    plt.close()
    print(f"  Saved: New_Human_Mouse_Comparison.pdf (r={r:.3f})")
else:
    print("  Skipped: Missing Human or Mouse data")

# ============== PLOT 3: CELL TYPE CORRELATION NETWORK ==============
# REMOVED per user request
print("\nPlot 3: Cell Type Correlation Matrix - SKIPPED")

# ============== PLOT 4: SAMPLE VARIABILITY RAINCLOUD ==============
print("\nPlot 4: Sample Variability Raincloud...")

# Select top 6 cell types by mean fraction
top_cts = combined_df.groupby('CellType')['Fraction'].mean().nlargest(6).index

fig, axes = plt.subplots(2, 3, figsize=(10, 7))
axes = axes.flatten()

for ax, ct in zip(axes, top_cts):
    ct_data = combined_df[combined_df['CellType'] == ct]
    
    # Create raincloud-like plot (violin + strip)
    # Map dataset names to descriptive labels
    ct_data = ct_data.copy()
    ct_data['DatasetLabel'] = ct_data['Dataset'].map(DATASETS)
    
    sns.violinplot(data=ct_data, x='DatasetLabel', y='Fraction', 
                  ax=ax, inner=None, alpha=0.6, cut=0,
                  palette=[COLORS['primary'], COLORS['secondary'], 
                          COLORS['tertiary'], COLORS['quaternary'], COLORS['background']])
    sns.stripplot(data=ct_data, x='DatasetLabel', y='Fraction', 
                 ax=ax, size=3, alpha=0.5, color='black')
    
    ax.set_title(ct, fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Fraction')
    ax.tick_params(axis='x', rotation=45)

plt.suptitle(f'Sample-Level Variability by Cell Type ({METHOD})', fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(OUTPUT_DIR / f"{METHOD}_New_Sample_Variability_Raincloud.pdf", bbox_inches='tight')
plt.close()
print("  Saved: New_Sample_Variability_Raincloud.pdf")

# ============== PLOT 5: RADAR CHART ==============
# REMOVED per user request
print("\nPlot 5: Dataset Comparison Radar Chart - SKIPPED")

# ============== SUMMARY ==============
print("\n" + "="*60)
print("Extended Deconvolution Analysis Complete!")
print("="*60)
print(f"\nNew plots saved to: {OUTPUT_DIR}")
print("\nGenerated:")
print("  1. New_Disease_Progression_Trajectory.pdf")
print("  2. New_Human_Mouse_Comparison.pdf")
# print("  3. New_CellType_Correlation_Matrix.pdf") (REMOVED)
print("  4. New_Sample_Variability_Raincloud.pdf")
# print("  5. New_Dataset_Radar_Chart.pdf") (REMOVED)
