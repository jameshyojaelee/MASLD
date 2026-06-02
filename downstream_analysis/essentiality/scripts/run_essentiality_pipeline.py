#!/usr/bin/env python3
"""
Master Essentiality Analysis Pipeline
======================================

This script runs the complete essentiality analysis pipeline, generating all plots
and CSV outputs. It is designed to be reproducible with any DEG list.

Usage:
    python run_essentiality_pipeline.py [--deg-file PATH]
    
    If --deg-file is not specified, uses updated_core_degs.csv by default.
    
Outputs:
    - 20+ publication-ready plots (PNG and PDF)
    - Quadrant gene lists (Essential_Both, Liang_Only, DepMap_Only, Non_Essential)
    - Combined essentiality CSVs

Author: Auto-generated for Sanjana Lab
"""

import argparse
import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Force non-interactive backend for HPC
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================================
# CONFIGURATION
# ============================================================================
BASE_DIR = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
ESSENTIALITY_DIR = BASE_DIR / "downstream_analysis/essentiality"
OUTPUT_DIR = ESSENTIALITY_DIR / "plots"
DEFAULT_DEG_FILE = BASE_DIR / "final_core_degs.csv"

# Reference data paths
LIANG_ESSENTIALITY_FILE = BASE_DIR / "RNA-seq/reference/extracted_tables/TableS2.xlsx"
ORTHOLOG_FILE = BASE_DIR / "streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"
MASTER_MATRIX_FILE = BASE_DIR / "streamlit_deg_explorer/master_ortholog_matrix.csv.gz"
DEPMAP_EFFECT_FILE = ESSENTIALITY_DIR / "CRISPRGeneEffect.csv"
DEPMAP_MODEL_FILE = ESSENTIALITY_DIR / "Model.csv"

# Liang screen sheets
LIANG_SHEETS = {
    "S2F": "HAP1",
    "S2G": "HEK293FT", 
    "S2H": "K562",
    "S2I": "MDA-MB-231",
    "S2J": "THP1"
}

# Thresholds
ESSENTIAL_THRESH = -0.5

# Cell lines to exclude from analysis
EXCLUDED_CELL_LINES = ['KMCH1', 'SKHEP1']

# Colors - Sanjana Lab palette
COLOR_CORE = "#e14b9d"  # Vivid Magenta/Pink
COLOR_BG = "#9c9c9c"    # Medium Gray

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def setup_plotting():
    """Configure matplotlib for publication-quality plots."""
    # Base context (paper)
    sns.set_context("paper") 
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
    
def format_ax(ax):
    """Standard axis formatting."""
    sns.despine(ax=ax, offset=10, trim=True)
    ax.grid(True, axis='y', linestyle='--', alpha=0.5, color='lightgray')
    ax.grid(True, axis='x', linestyle='--', alpha=0.5, color='lightgray')
    ax.tick_params(axis='both', length=6)

def get_origin(analysis_str):
    """Classify gene origin from analysis string."""
    a = str(analysis_str).lower()
    has_mouse = 'mouse' in a or 'mcd' in a
    has_human = 'patient' in a or 'gwas' in a or 'human' in a
    if has_mouse and has_human:
        return 'Both'
    elif has_mouse:
        return 'Mouse'
    elif has_human:
        return 'Human'
    return 'Other'

def get_quadrant(row, liang_col='Mean_Essentiality_Score', depmap_col='DepMap_Cas9_Score'):
    """Classify gene into essentiality quadrant."""
    liang_ess = row[liang_col] < ESSENTIAL_THRESH
    depmap_ess = row[depmap_col] < ESSENTIAL_THRESH
    if liang_ess and depmap_ess:
        return 'Essential_Both'
    elif liang_ess:
        return 'Essential_Liang_Only'
    elif depmap_ess:
        return 'Essential_DepMap_Only'
    else:
        return 'Non_Essential'

# ============================================================================
# DATA LOADING
# ============================================================================
def load_core_degs(deg_file):
    """
    Load Core DEG file and adapt Mouse-Centric format to Human-Centric analysis.
    The script expects 'gene_symbol', 'species', 'analyses', 'overlap_count'.
    """
    print(f"Loading Core DEGs from {deg_file}...")
    df = pd.read_csv(deg_file)
    print(f"  Loaded {len(df)} mouse targets")
    
    # ADAPTER: final_core_degs.csv (Mouse Centric) -> Analysis Format (Human)
    # 1. We need Human Symbols for DepMap/Liang checking.
    # 2. Key columns: human_ortholog_symbols, source_analyses
    
    # Filter for targets that HAVE a human ortholog
    # (Mouse-specific targets cannot be validated in Human DepMap)
    df_human = df.dropna(subset=['human_ortholog_symbols']).copy()
    df_human = df_human[df_human['human_ortholog_symbols'] != ""]
    
    print(f"  Targets with Human Orthologs: {len(df_human)}")
    
    # Explode if multiple orthologs (separated by ;)
    # We want to check ALL human orthologs associated with our library
    expanded_rows = []
    
    for _, row in df_human.iterrows():
        h_syms = str(row['human_ortholog_symbols']).split(';')
        h_ids = str(row['human_ortholog_ids']).split(';')
        
        # Ensure lengths match or handle gracefully
        # Usually they should match.
        
        for hs in h_syms:
            if hs.strip():
                expanded_rows.append({
                    'gene_symbol': hs.strip(),
                    'mouse_gene_id': row['mouse_gene_id'], # Keep track of source
                    'analyses': row['source_analyses'],
                    'overlap_count': row['n_analyses'], # Use n_analyses as proxy for overlap
                    'species': 'human' # Mock species to pass downstream filters
                })
                
    df_adapted = pd.DataFrame(expanded_rows)
    print(f"  Expanded Human Orthologs: {len(df_adapted)}")
    
    if len(df_adapted) > 0:
        print(f"  Unique Human Genes: {df_adapted['gene_symbol'].nunique()}")
        
    return df_adapted

def load_liang_data():
    """Load Liang et al. Cas13 screen data."""
    print(f"Loading Liang Cas13 screen data...")
    screen_data = pd.DataFrame()
    
    for sheet, cell_line in LIANG_SHEETS.items():
        df = pd.read_excel(LIANG_ESSENTIALITY_FILE, sheet_name=sheet, header=2, engine='openpyxl')
        df.columns = [str(c).strip() for c in df.columns]
        
        # Find gene and FC columns
        gene_col = next((c for c in df.columns if 'Gene' in c or 'id' in c.lower()), None)
        fc_cols = [c for c in df.columns if 'fold-change' in c.lower() or 'fc' in c.lower()]
        
        if gene_col and fc_cols:
            day14_col = next((c for c in fc_cols if '14' in c), fc_cols[0])
            subset = df[[gene_col, day14_col]].copy()
            subset.columns = ['Gene', 'Score']
            subset['CellLine'] = cell_line
            screen_data = pd.concat([screen_data, subset])
    
    print(f"  Loaded {len(screen_data)} scores from {len(LIANG_SHEETS)} cell lines")
    return screen_data

def load_depmap_liver_data():
    """Load DepMap Cas9 screen data for liver cell lines."""
    print(f"Loading DepMap Cas9 screen data...")
    
    model_df = pd.read_csv(DEPMAP_MODEL_FILE)
    lineage_col = 'OncotreeLineage' if 'OncotreeLineage' in model_df.columns else None
    
    if not lineage_col:
        print("  ERROR: Could not find lineage column in Model.csv")
        return None, None
    
    # Get liver models
    liver_models_df = model_df[model_df[lineage_col] == 'Liver']
    
    # Exclude specified cell lines
    name_col = 'StrippedCellLineName' if 'StrippedCellLineName' in model_df.columns else 'CellLineName'
    excluded_mask = liver_models_df[name_col].isin(EXCLUDED_CELL_LINES)
    n_excluded = excluded_mask.sum()
    liver_models_df = liver_models_df[~excluded_mask]
    print(f"  Excluded {n_excluded} cell lines: {EXCLUDED_CELL_LINES}")
    
    liver_models = liver_models_df['ModelID'].tolist()
    id_to_name = liver_models_df.set_index('ModelID')[name_col].to_dict()
    
    print(f"  Found {len(liver_models)} liver cell lines (after exclusion)")
    
    # Load effect data
    full_effect = pd.read_csv(DEPMAP_EFFECT_FILE, index_col=0)
    liver_effect = full_effect[full_effect.index.isin(liver_models)]
    liver_effect.rename(index=id_to_name, inplace=True)
    liver_effect = liver_effect.groupby(level=0).mean()
    
    # Clean column names
    liver_effect.columns = [c.split(' ')[0] for c in liver_effect.columns]
    liver_effect = liver_effect.loc[:, ~liver_effect.columns.duplicated()]
    
    return liver_effect, len(liver_models)

# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================
def run_pipeline(deg_file):
    """Run the complete essentiality analysis pipeline."""
    
    setup_plotting()
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # -------------------------------------------------------------------------
    # 1. LOAD DATA
    # -------------------------------------------------------------------------
    core_df = load_core_degs(deg_file)
    
    # Get human gene symbols (for matching with screens)
    # The adapter sets 'species'='human' for all rows
    human_rows = core_df[core_df['species'] == 'human']
    foreground_symbols = set(human_rows['gene_symbol'].unique())
    print(f"  Foreground contains {len(foreground_symbols)} unique Human symbols")
    
    # Add metadata
    # We aggregate if multiple mouse genes map to same human gene
    # Keep max overlap count to be generous
    meta_df = core_df[['gene_symbol', 'overlap_count', 'analyses']].copy()
    meta_df = meta_df.sort_values('overlap_count', ascending=False).drop_duplicates(subset='gene_symbol')
    
    meta_df['Origin'] = meta_df['analyses'].apply(get_origin)
    
    # Load Liang data
    liang_data = load_liang_data()
    mean_liang_scores = liang_data.groupby('Gene')['Score'].mean().reset_index()
    mean_liang_scores.columns = ['Gene', 'Mean_Essentiality_Score']
    
    # Separate foreground and background
    foreground_liang = mean_liang_scores[mean_liang_scores['Gene'].isin(foreground_symbols)]
    background_liang = mean_liang_scores[~mean_liang_scores['Gene'].isin(foreground_symbols)]
    
    print(f"Foreground genes in Liang screen: {len(foreground_liang)}")
    print(f"Background genes in Liang screen: {len(background_liang)}")
    
    # Merge metadata
    merged_foreground = foreground_liang.merge(meta_df, left_on='Gene', right_on='gene_symbol', how='left')
    
    # Load DepMap data
    liver_effect, n_liver_models = load_depmap_liver_data()
    
    if liver_effect is not None:
        # Calculate DepMap metrics
        target_genes = set(merged_foreground['Gene'])
        available_genes = list(set(liver_effect.columns) & target_genes)
        print(f"Core genes in DepMap: {len(available_genes)}")
        
        liver_effect_core = liver_effect[available_genes]
        depmap_mean = liver_effect_core.mean(axis=0)
        depmap_mean.name = 'DepMap_Cas9_Score'
        depmap_overlap = (liver_effect_core < ESSENTIAL_THRESH).sum(axis=0)
        depmap_overlap.name = 'DepMap_Overlap_Count'
        
        merged_foreground = merged_foreground.merge(depmap_mean, left_on='Gene', right_index=True, how='left')
        merged_foreground = merged_foreground.merge(depmap_overlap, left_on='Gene', right_index=True, how='left')
        
        # Calculate percentages
        merged_foreground['Liang_Pct_Essential'] = merged_foreground['overlap_count'] / 5 * 100
        merged_foreground['DepMap_Pct_Essential'] = merged_foreground['DepMap_Overlap_Count'] / n_liver_models * 100
    
    # -------------------------------------------------------------------------
    # 2. GENERATE PLOTS
    # -------------------------------------------------------------------------
    print("\n--- Generating Plots ---")
    
    # -------------------------------------------------------------------------
    # DepMap Unified Plots (if available)
    # -------------------------------------------------------------------------
    if liver_effect is not None and 'DepMap_Cas9_Score' in merged_foreground.columns:
        print("  3. Unified Scatter Comparison (Advanced)...")
        # 1. Prepare Main Data (Quadrants)
        plot_df = merged_foreground.dropna(subset=['Mean_Essentiality_Score', 'DepMap_Cas9_Score']).copy()
        plot_df['Quadrant'] = plot_df.apply(get_quadrant, axis=1)
        
        # 2. Prepare Side Panel Data (Not in Liang)
        screened_genes = set(plot_df['Gene'])
        missing_genes = list(foreground_symbols - screened_genes)
        
        side_panel_df = pd.DataFrame()
        if missing_genes:
            # Get DepMap scores for these genes
            # Use full liver_effect to check for genes NOT in core subset
            valid_missing = [g for g in missing_genes if g in liver_effect.columns]
            if valid_missing:
                missing_scores = liver_effect[valid_missing].mean(axis=0)
                side_panel_df = pd.DataFrame({'DepMap_Cas9_Score': missing_scores, 'Gene': missing_scores.index})
                # Empty string to avoid giant tick label
                side_panel_df['Group'] = ''
        
        # 3. Plotting
        # Custom Palette
        palette = {
            'Non_Essential': '#e0e0e0',          # Light Grey
            'Essential_Liang_Only': '#f3c347',   # Yellow/Gold
            'Essential_DepMap_Only': '#ff7f0e',  # Orange
            'Essential_Both': '#e14b9d'          # Pink (Core)
        }
        
        # Sort for z-order (Non-essential first)
        plot_df['Z'] = plot_df['Quadrant'].map({
            'Non_Essential': 0, 
            'Essential_Liang_Only': 1, 
            'Essential_DepMap_Only': 1, 
            'Essential_Both': 2
        })
        plot_df = plot_df.sort_values('Z')
        
        # Setup Figure
        fig, (ax_main, ax_side) = plt.subplots(1, 2, figsize=(8, 5), width_ratios=[5, 1], sharey=True)
        plt.subplots_adjust(wspace=0.1)
        # Main Scatter
        # Prepare Legend Labels (Simplify and Add Counts)
        counts = plot_df['Quadrant'].value_counts()
        name_map = {
            'Essential_Both': 'Both',
            'Essential_Liang_Only': 'Liang Only',
            'Essential_DepMap_Only': 'DepMap Only',
            'Non_Essential': 'Non-Essential'
        }
        
        # Create mapping: Internal -> "Display Name (n=X)"
        legend_map = {}
        for q in name_map: # Ensure we cover all keys even if count is 0
            count = counts.get(q, 0)
            display = name_map[q]
            legend_map[q] = f"{display} (n={count})"
            
        plot_df['Legend_Label'] = plot_df['Quadrant'].map(legend_map)

        # Update Palette to match new Legend Labels
        # Original Palette logic:
        # palette = {'Essential_Both': '#e14b9d', ...}
        # We need to map new keys to these values
        base_colors = {
            'Both': '#e14b9d',
            'Liang Only': '#FBC02D',
            'DepMap Only': '#FF7043', 
            'Non-Essential': '#E0E0E0'
        }
        new_palette = {}
        for internal, new_label in legend_map.items():
            simple_name = name_map.get(internal)
            if simple_name in base_colors:
                new_palette[new_label] = base_colors[simple_name]

        # Main Scatter
        sns.scatterplot(data=plot_df, x='Mean_Essentiality_Score', y='DepMap_Cas9_Score',
                        hue='Legend_Label', palette=new_palette, 
                        alpha=0.8, s=30, edgecolor='white', linewidth=0.3, ax=ax_main)
        
        # Threshold Lines
        ax_main.axhline(ESSENTIAL_THRESH, linestyle='--', color='lightgray', linewidth=1)
        ax_main.axvline(ESSENTIAL_THRESH, linestyle='--', color='lightgray', linewidth=1)
        
        ax_main.set_title('Cross-Modality Essentiality Comparison', fontsize=16, fontweight='bold')
        ax_main.set_xlabel('Liang Cas13 Score (Mean Log2FC)', fontsize=14)
        ax_main.set_ylabel('DepMap Cas9 Score (Mean)', fontsize=14)
        # Legend (Enlarged)
        ax_main.legend(title='Shared Genes Status', loc='upper left', frameon=True,
                       fontsize=12, title_fontsize=13, framealpha=0.95, borderpad=1)
        
        # Removed manual text annotations on plot as requested
                         
        format_ax(ax_main)
        
        # Side Panel (Strip Plot)
        if not side_panel_df.empty:
            # Color points based on essentiality in DepMap
            side_panel_df['Is_Essential'] = side_panel_df['DepMap_Cas9_Score'] < ESSENTIAL_THRESH
            side_colors = {True: '#ff7f0e', False: '#e0e0e0'} # Orange if essential, Grey if not
            
            # We use stripplot for jitter
            sns.stripplot(data=side_panel_df, y='DepMap_Cas9_Score', x='Group',
                          hue='Is_Essential', palette=side_colors,
                          alpha=0.6, s=4, jitter=0.25, ax=ax_side, legend=False)
            
            # Threshold Line
            ax_side.axhline(ESSENTIAL_THRESH, linestyle='--', color='lightgray', linewidth=1)
            
            # Annotate counts
            n_missing = len(side_panel_df)
            ax_side.set_xlabel(f'Not Screened\nin Liang\n(n={n_missing})', fontsize=12)
            # Annotation
            n_ess = side_panel_df['Is_Essential'].sum()
            if n_ess > 0:
                ax_side.text(0, min(side_panel_df['DepMap_Cas9_Score']) - 0.2, 
                             f'n={n_ess}', ha='center', va='top', 
                             fontweight='bold', fontsize=10, 
                             bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=2))
        
        ax_side.set_ylabel('') # Hide Y label
        format_ax(ax_side)
        ax_side.grid(False, axis='x') # Keep only Y grid
        
        # Save
        plt.tight_layout()
        plt.savefig(OUTPUT_DIR / 'Unified_Scatter_Comparison.pdf', dpi=300, bbox_inches='tight')
        plt.close()
    
    # -------------------------------------------------------------------------
    # 3. GENERATE QUADRANT EXPORTS
    # -------------------------------------------------------------------------
    print("\n--- Exporting Quadrant Gene Lists ---")
    
    RESULTS_DIR = ESSENTIALITY_DIR / "results"
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    if 'DepMap_Cas9_Score' in merged_foreground.columns:
        df_valid = merged_foreground.dropna(subset=['Mean_Essentiality_Score', 'DepMap_Cas9_Score']).copy()
        df_valid['Essentiality_Quadrant'] = df_valid.apply(get_quadrant, axis=1)
        
        # Standard Quadrants
        for quadrant in ['Essential_Both', 'Essential_Liang_Only', 'Essential_DepMap_Only', 'Non_Essential']:
            subset = df_valid[df_valid['Essentiality_Quadrant'] == quadrant].copy()
            subset = subset.sort_values('Mean_Essentiality_Score')
            cols = ['Gene', 'Mean_Essentiality_Score', 'DepMap_Cas9_Score', 'Origin', 'overlap_count', 'analyses']
            subset[cols].to_csv(RESULTS_DIR / f'Quadrant_{quadrant}.csv', index=False)
            print(f"  {quadrant}: {len(subset)} genes")
        
        # Requested Specific Files
        print("\n--- Generating Requested User Files ---")
        
        # 1. Essential in Both
        essential_both = df_valid[df_valid['Essentiality_Quadrant'] == 'Essential_Both'].sort_values('Mean_Essentiality_Score')
        essential_both[cols].to_csv(RESULTS_DIR / '1_Essential_Both_DepMap_Liang.csv', index=False)
        print(f"  1. Essential Both: {len(essential_both)}")
        
        # 2. Essential in DepMap (All) -> Union of Both + DepMap_Only
        # Logic: DepMap_Cas9_Score < THRESH (-0.5)
        essential_depmap_all = df_valid[df_valid['DepMap_Cas9_Score'] < ESSENTIAL_THRESH].sort_values('DepMap_Cas9_Score')
        essential_depmap_all[cols].to_csv(RESULTS_DIR / '2_Essential_DepMap_All.csv', index=False)
        print(f"  2. Essential DepMap (All): {len(essential_depmap_all)}")
        
        # 3. Essential in Liang Only
        essential_liang_only = df_valid[df_valid['Essentiality_Quadrant'] == 'Essential_Liang_Only'].sort_values('Mean_Essentiality_Score')
        essential_liang_only[cols].to_csv(RESULTS_DIR / '3_Essential_Liang_Only.csv', index=False)
        print(f"  3. Essential Liang Only: {len(essential_liang_only)}")
        
        # Combined file
        df_valid_export = df_valid[['Gene', 'Mean_Essentiality_Score', 'DepMap_Cas9_Score', 
                                     'Origin', 'overlap_count', 'Essentiality_Quadrant', 'analyses']].copy()
        df_valid_export.to_csv(RESULTS_DIR / 'Essentiality_Quadrants_All.csv', index=False)
        print(f"  Combined: {len(df_valid_export)} genes -> Essentiality_Quadrants_All.csv")
    
    # -------------------------------------------------------------------------
    # 4. SAVE MAIN CSV OUTPUTS
    # -------------------------------------------------------------------------
    print("\n--- Saving CSV Outputs ---")
    foreground_liang.to_csv(RESULTS_DIR / 'Core_DEGs_Essentiality.csv', index=False)
    print(f"  Core_DEGs_Essentiality.csv: {len(foreground_liang)} genes")
    
    merged_foreground.to_csv(RESULTS_DIR / 'Unified_Essentiality_Combined.csv', index=False)
    print(f"  Unified_Essentiality_Combined.csv: {len(merged_foreground)} genes")
    
    # -------------------------------------------------------------------------
    # 5. RUN EXTENDED ANALYSIS
    # -------------------------------------------------------------------------
    print("\n--- Running Extended Analysis ---")
    run_extended_analysis(merged_foreground, liver_effect, liver_effect_core if liver_effect is not None else None)
    
    print("\n=== Pipeline Complete ===")
    print(f"Outputs saved to: {OUTPUT_DIR}")


def run_extended_analysis(merged_df, liver_effect, liver_effect_core):
    """Run the extended analysis plots (Liver Heatmap)."""
    
    if liver_effect is None:
        print("  Skipping extended analysis (no DepMap data)")
        return
    
    # Mean essentiality per gene (all genes in DepMap)
    mean_essentiality = liver_effect.mean(axis=0).reset_index()
    mean_essentiality.columns = ['Gene', 'Mean_Essentiality']
    mean_essentiality = mean_essentiality.sort_values('Mean_Essentiality')
    mean_essentiality['Rank'] = range(1, len(mean_essentiality) + 1)
    
    core_symbols = set(merged_df['Gene'].dropna())
    core_in_depmap = mean_essentiality[mean_essentiality['Gene'].isin(core_symbols)]
    essential_core = core_in_depmap[core_in_depmap['Mean_Essentiality'] < ESSENTIAL_THRESH].copy()
    essential_core = essential_core.sort_values('Mean_Essentiality').head(50)
    
    # PLOT: Liver-Specific Heatmap
    print("  6. Liver Essentiality Heatmap...")
    if liver_effect_core is not None and len(liver_effect_core.columns) > 0:
        top_core_genes = essential_core.head(30)['Gene'].tolist() if len(essential_core) > 0 else []
        
        if top_core_genes:
            available_top = [g for g in top_core_genes if g in liver_effect_core.columns]
            if available_top:
                liver_heatmap = liver_effect_core[available_top].T
                liver_heatmap['mean'] = liver_heatmap.mean(axis=1)
                liver_heatmap = liver_heatmap.sort_values('mean').drop('mean', axis=1)
                
                plt.figure(figsize=(8, 6))
                sns.heatmap(liver_heatmap, cmap='RdBu', center=0, 
                            cbar_kws={'label': 'Essentiality Score'},
                            linewidths=0.5)
                plt.title('Core DEG Essentiality in Liver Cell Lines (DepMap)', fontsize=10, fontweight='bold')
                plt.xlabel('Liver Cell Line', fontsize=8)
                plt.ylabel('Gene', fontsize=8)
                plt.tight_layout()
                # PNG removed
                plt.savefig(OUTPUT_DIR / 'New_Liver_Essentiality_Heatmap.pdf', bbox_inches='tight')
                plt.close()


# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run essentiality analysis pipeline')
    parser.add_argument('--deg-file', type=str, default=str(DEFAULT_DEG_FILE),
                        help=f'Path to DEG file (default: {DEFAULT_DEG_FILE})')
    args = parser.parse_args()
    
    deg_file = Path(args.deg_file)
    if not deg_file.exists():
        print(f"ERROR: DEG file not found: {deg_file}")
        sys.exit(1)
    
    run_pipeline(deg_file)
