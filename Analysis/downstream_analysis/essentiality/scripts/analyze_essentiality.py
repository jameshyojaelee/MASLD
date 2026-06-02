import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys

# --- Configuration ---
CORE_DEGS_FILE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/final_core_degs.csv"
ESSENTIALITY_FILE = "RNA-seq/reference/extracted_tables/TableS2.xlsx"
ORTHOLOG_FILE = "streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"
OUTPUT_DIR = "Analysis/downstream_analysis/essentiality/plots/"
SHEETS_TO_PROCESS = {
    "S2F": "HAP1",
    "S2G": "HEK293FT",
    "S2H": "K562",
    "S2I": "MDA-MB-231",
    "S2J": "THP1"
}

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Starting Essentiality Analysis...")

# --- 1. Load Core DEGs ---
# --- 1. Load Core DEGs (with Adapter for Mouse-Centric Format) ---
print(f"Loading Core DEGs from {CORE_DEGS_FILE}...")
raw_df = pd.read_csv(CORE_DEGS_FILE)
print(f"  Raw shape: {raw_df.shape}")

# ADAPTER: final_core_degs.csv -> Standard Analysis Format
# Filter for targets that HAVE a human ortholog
df_human = raw_df.dropna(subset=['human_ortholog_symbols']).copy()
df_human = df_human[df_human['human_ortholog_symbols'] != ""]

rows = []
for _, row in df_human.iterrows():
    h_syms = str(row['human_ortholog_symbols']).split(';')
    for hs in h_syms:
        if hs.strip():
            rows.append({
                'gene_symbol': hs.strip(),
                'mouse_gene_id': row['mouse_gene_id'],
                'analyses': row['source_analyses'],
                'overlap_count': row['n_analyses'],
                'species': 'human', # Treat as human for essentiality lookup
                'human_ortholog_gene': hs.strip() # Mock for backward compat
            })

core_df = pd.DataFrame(rows)
print(f"  Adapted to {len(core_df)} human ortholog rows")
print(f"  Unique Human Genes: {core_df['gene_symbol'].nunique()}")

# Use these mapped IDs as our "Foreground".
foreground_genes = set(core_df['gene_symbol'].unique())
mapped_human_ids = foreground_genes # Alias for legacy code usage if any

# Use these mapped IDs as our "Foreground".
foreground_genes = set(mapped_human_ids)

# --- 3. Load Essentiality Data ---
print(f"Loading Essentiality Data from {ESSENTIALITY_FILE}...")
essentiality_scores = []

# Collect scores from all cell lines
full_screen_genes = set()

# Define Foreground Symbols
# The user's input CSV has `human` rows too!
human_rows = core_df[core_df['species'] == 'human']
print(f"Found {len(human_rows)} Human rows in input CSV.")
foreground_symbols = set(human_rows['gene_symbol'].unique())
print(f"Foreground contains {len(foreground_symbols)} unique Human symbols from input CSV.")

# If we have 0 human rows, we must do the mapping.
if len(human_rows) == 0:
    print("WARNING: No human rows found. We need Ensembl -> Symbol mapping.")
    pass

# Load Screen Data
screen_data = pd.DataFrame()

for sheet, cell_line in SHEETS_TO_PROCESS.items():
    print(f"  Processing {sheet} ({cell_line})...")
    # Read sheet. Skip header rows if needed.
    # header=2 means the 3rd row is the header (Index 2).
    # Based on inspection: Row 0=Title, Row 1=NaNs, Row 2=Header ('Gene', etc)
    df = pd.read_excel(ESSENTIALITY_FILE, sheet_name=sheet, header=2, engine='openpyxl')
    
    # Check columns
    # We want 'Gene' (Symbol) or ID?
    df.columns = [str(c).strip() for c in df.columns]
    
    # Find Gene column
    gene_col = [c for c in df.columns if 'Gene' in c or 'id' in c.lower()]
    if not gene_col:
        print(f"    Could not find Gene column in {sheet}. Cols: {df.columns}")
        # Fallback: check if 'Gene' is in the first row of data if header was wrong?
        # But we trust header=2 for now.
        continue
    gene_col = gene_col[0]

    
    # Find FC column
    # Often 'Fold-change' or 'log2FC'
    fc_cols = [c for c in df.columns if 'fold-change' in c.lower() or 'fc' in c.lower()]
    if not fc_cols:
         print(f"    Could not find FC column in {sheet}. Cols: {df.columns}")
         continue
    # Prefer Day 14 or latest if multiple?
    # Paper usually uses Day 14 or Average. Let's take the first one found or specifically Day 14 if available.
    day14_col = [c for c in fc_cols if '14' in c]
    target_fc_col = day14_col[0] if day14_col else fc_cols[0]
    
    print(f"    Using Gene Col: {gene_col}, Score Col: {target_fc_col}")
    
    subset = df[[gene_col, target_fc_col]].copy()
    subset.columns = ['Gene', 'Score']
    subset['CellLine'] = cell_line
    screen_data = pd.concat([screen_data, subset])

print(f"Loaded {len(screen_data)} scores.")

# --- 4. Analysis ---
# Calculate Mean Essentiality across cell lines for each gene
mean_scores = screen_data.groupby('Gene')['Score'].mean().reset_index()
mean_scores.columns = ['Gene', 'Mean_Essentiality_Score']

# Define Sets
# Foreground: Genes in our input list (Symbols)
# Background: All other genes in the screen
# Intersect foreground with screen genes
foreground_in_screen = mean_scores[mean_scores['Gene'].isin(foreground_symbols)]
background_in_screen = mean_scores[~mean_scores['Gene'].isin(foreground_symbols)]

print(f"Foreground genes in screen: {len(foreground_in_screen)}")
print(f"Background genes in screen: {len(background_in_screen)}")

if len(foreground_in_screen) == 0:
    print("ERROR: No overlap between Input Symbols and Screen Symbols.")
    print("Input Sample:", list(foreground_symbols)[:5])
    print("Screen Sample:", mean_scores['Gene'].head().tolist())
    # Try mapping Mouse -> Human if input was mouse-only?
    # But for now assume input has human symbols.
else:
    # --- 5. Plotting ---
    print("Generating Plots...")
    
    # Merge Metadata for Creative Plots
    # We need 'overlap_count' and 'analyses' from core_df mapped to foreground_in_screen
    meta_df = core_df[['gene_symbol', 'overlap_count', 'analyses']].drop_duplicates(subset='gene_symbol')
    
    # Define Origin Helper
    def get_origin(analysis_str):
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

    meta_df['Origin'] = meta_df['analyses'].apply(get_origin)
    
    merged_foreground = foreground_in_screen.merge(meta_df, left_on='Gene', right_on='gene_symbol', how='left')
    
    # --- 5. Plotting (Enhanced & Creative) ---
    print("Generating Comprehensive Plots...")
    
    # --- Style Configuration ---
    COLOR_CORE = "#e14b9d"  # Vivid Magenta/Pink
    COLOR_BG = "#9c9c9c"    # Medium Gray
    COLOR_PALETTE = [COLOR_CORE, COLOR_BG]
    
    sns.set_context("paper", font_scale=2.0)
    sns.set_style("whitegrid", {'axes.grid': False})
    
    # Configure fonts for Illustrator (Type 42 = TrueType)
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # Font Family (Helvetica/Arial)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'sans-serif']

    # Strict Font Sizes (Journal Standards)
    plt.rcParams['font.size'] = 7           # Default text size
    plt.rcParams['axes.titlesize'] = 8      # Title size
    plt.rcParams['axes.labelsize'] = 8      # X and Y label size
    plt.rcParams['xtick.labelsize'] = 6     # X tick size
    plt.rcParams['ytick.labelsize'] = 6     # Y tick size
    plt.rcParams['legend.fontsize'] = 6     # Legend size
    
    def format_ax(ax):
        sns.despine(ax=ax, offset=10, trim=True)
        ax.grid(True, axis='y', linestyle='--', alpha=0.5, color='lightgray')
        ax.grid(True, axis='x', linestyle='--', alpha=0.5, color='lightgray')
        ax.tick_params(axis='both', length=6)

    # ... [CDF, Density, Raincloud code remains same or similar] ...
    # 1. CDF Plot - REMOVED per user request
    # plt.figure(figsize=(10, 7))
    # ax = plt.gca()
    # sns.ecdfplot(data=background_in_screen, x='Mean_Essentiality_Score', 
    #              label=f'Background (n={len(background_in_screen)})', 
    #              color=COLOR_BG, linewidth=2, alpha=0.8)
    # sns.ecdfplot(data=foreground_in_screen, x='Mean_Essentiality_Score', 
    #              label=f'Core PCG DEGs (n={len(foreground_in_screen)})', 
    #              color=COLOR_CORE, linewidth=3)
    # plt.title('Gene Essentiality CDF', fontsize=18, fontweight='bold', pad=20)
    # plt.xlabel('Essentiality Score (Mean Log2FC)', fontsize=14)
    # plt.ylabel('Cumulative Probability', fontsize=14)
    # plt.legend(frameon=False, fontsize=12)
    # format_ax(ax)
    # plt.tight_layout()
    # plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_CDF.pdf'), dpi=300)
    # plt.close()
    
    # 2. Density Plot
    plt.figure(figsize=(10, 7))
    ax = plt.gca()
    sns.kdeplot(data=background_in_screen, x='Mean_Essentiality_Score', 
                label='Background', color=COLOR_BG, fill=True, alpha=0.2, linewidth=0)
    sns.kdeplot(data=foreground_in_screen, x='Mean_Essentiality_Score', 
                label='Core PCG DEGs', color=COLOR_CORE, fill=True, alpha=0.4, linewidth=2)
    plt.title('Essentiality Density Profile', fontsize=18, fontweight='bold', pad=20)
    plt.xlabel('Essentiality Score (Mean Log2FC)', fontsize=14)
    plt.ylabel('Density', fontsize=14)
    plt.legend(frameon=False, fontsize=12)
    format_ax(ax)
    sns.despine(left=True)
    ax.yaxis.set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Density.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Density.pdf'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Raincloud Plot
    df_core = pd.DataFrame({'Score': foreground_in_screen['Mean_Essentiality_Score'], 'Group': 'Core PCG DEGs'})
    df_bg = pd.DataFrame({'Score': background_in_screen['Mean_Essentiality_Score'], 'Group': 'Background'})
    if len(df_bg) > 1000:
        df_bg_strip = df_bg.sample(1000, random_state=42)
    else:
        df_bg_strip = df_bg
    
    plot_data = pd.concat([df_core, df_bg])
    plot_data_strip = pd.concat([df_core, df_bg_strip])
    
    plt.figure(figsize=(8, 8))
    ax = plt.gca()
    sns.boxplot(data=plot_data, x='Group', y='Score', 
                palette=[COLOR_CORE, COLOR_BG], 
                width=0.4, fliersize=0, linewidth=1.5,
                boxprops=dict(alpha=0.7))
    sns.stripplot(data=plot_data_strip, x='Group', y='Score', 
                  palette=[COLOR_CORE, COLOR_BG], 
                  size=2, alpha=0.4, jitter=0.2)
    plt.title('Essentiality Distribution', fontsize=18, fontweight='bold', pad=20)
    plt.xlabel('')
    plt.ylabel('Essentiality Score (Mean Log2FC)', fontsize=14)
    format_ax(ax)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Raincloud.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Raincloud.pdf'), dpi=300, bbox_inches='tight')
    plt.close()

    # 4. OVERVIEW SUNBURST (Implemented as Nested Pie due to library constraints)
    # Layer 1: Origin (Mouse, Human, Both)
    # Layer 2: In Screen vs Not In Screen
    # Layer 3: Essential (<-0.5) vs Non-Essential
    
    # Prepare Data
    # All Core Genes (Total)
    all_core_meta = meta_df[['gene_symbol', 'Origin']].copy()
    all_core_meta['InScreen'] = all_core_meta['gene_symbol'].isin(foreground_in_screen['Gene'])
    
    # Merge Score to get Essentiality
    all_core_meta = all_core_meta.merge(foreground_in_screen[['Gene', 'Mean_Essentiality_Score']], 
                                        left_on='gene_symbol', right_on='Gene', how='left')
    
    # Define Essentiality Threshold (e.g., < -0.5 is Essential)
    ESSENTIAL_THRESH = -0.5
    all_core_meta['Status'] = 'Not screened'
    all_core_meta.loc[all_core_meta['InScreen'], 'Status'] = 'Non-Essential'
    all_core_meta.loc[(all_core_meta['InScreen']) & (all_core_meta['Mean_Essentiality_Score'] < ESSENTIAL_THRESH), 'Status'] = 'Essential'

    # Aggregations
    # L1: Origin
    l1_counts = all_core_meta.groupby('Origin').size()
    # L2: Origin + InScreen
    l2_counts = all_core_meta.groupby(['Origin', 'InScreen']).size()
    # L3: Origin + InScreen + Status
    l3_counts = all_core_meta.groupby(['Origin', 'InScreen', 'Status']).size()
    
    # Plotting Nested Donut
    plt.figure(figsize=(10, 10))
    ax = plt.gca()
    
    # Colors
    cmap = plt.get_cmap("tab20c")
    outer_colors = cmap(np.arange(3)*4)
    # Custom mapping might be better but let's try auto first
    
    # Ring 1: Origin
    pie_1_data = l1_counts
    plt.pie(pie_1_data, labels=pie_1_data.index, radius=1, 
            wedgeprops=dict(width=0.3, edgecolor='w'), pctdistance=0.85, 
            autopct='%1.1f%%', colors=['#ff9999','#66b3ff','#99ff99','#ffcc99'])
            
    # Ring 2: Coverage (Inside) - Logic roughly
    # Actually standard nested pie is hard to align perfectly without hierarchical data preparation
    # Let's simplify: Just Species breakdown of ESSENTIAL genes?
    # User asked for "Overview of how many genes in our core file are essential, distribution, coming from human/mouse"
    
    # Improved Sunburst alternative: Stacked Formatting
    # Let's stick to a cleaned up Nested Pie
    # Outer Ring: Origin
    # Inner Ring: Essentiality Status
    
    group_counts = all_core_meta.groupby(['Origin', 'Status']).size().reset_index(name='count')
    # Sort for alignment
    group_counts = group_counts.sort_values(['Origin', 'Status'])
    
    # Need to match colors. 
    # Let's do a simple 2-level pie:
    # Outer: Origin
    # Inner: Status
    
    # Re-calc for simple plotting
    # Pivot for stacked bar might be clearer? 
    # User specifically asked for creative comprehensive plots.
    # Let's do a "Hierarchical Bar" or "Mosaic Plot" logic
    # Or just 3 subplots?
    
    # Let's try the Nested Pie:
    # Outer: Origin
    # Inner: Essential vs Non-Essential (Screened only?)
    # Maybe just plot Screened genes?
    # No, user asked "comprehensive overview of how many genes in our core file..." implies ALL 3872.
    
    plt.clf()
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Stacked Bar with Data Table
    pivot_data = all_core_meta.groupby(['Origin', 'Status']).size().unstack(fill_value=0)
    # Reorder columns
    cols = ['Essential', 'Non-Essential', 'Not screened']
    # Ensure all columns exist
    for c in cols:
        if c not in pivot_data.columns:
            pivot_data[c] = 0
    pivot_data = pivot_data[cols]
    
    # Colors
    colors = ['#e14b9d', '#9c9c9c', '#e6e6e5'] # Pink, Gray, LightGray
    
    # 4. (Previously Stacked Bar - Removed in favor of Unified Side-by-Side version)

    
    # 5. Species-Split Violin Plot
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    
    sns.violinplot(data=merged_foreground, x='Origin', y='Mean_Essentiality_Score', 
                   hue='Origin', palette='viridis', cut=0, inner='quartile')
    
    plt.title('Essentiality Score Distribution by Origin', fontsize=18, fontweight='bold', pad=20)
    plt.xlabel('Gene Origin', fontsize=14)
    plt.ylabel('Essentiality Score (Mean Log2FC)', fontsize=14)
    format_ax(ax)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Species_Split_Violin.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Species_Split_Violin.pdf'), dpi=300, bbox_inches='tight')
    plt.close()


    # 5. Overlap Correlation Boxplot
    # ... [Overlap Boxplot code] ...
    # Ensure categorical order
    merged_foreground['overlap_group'] = merged_foreground['overlap_count'].apply(lambda x: str(int(x)) if pd.notnull(x) else 'Unknown')
    groups = sorted(merged_foreground['overlap_group'].unique())
    plt.figure(figsize=(8, 8))
    ax = plt.gca()
    sns.boxplot(data=merged_foreground, x='overlap_group', y='Mean_Essentiality_Score', 
                palette='magma_r', order=groups, width=0.5)
    sns.stripplot(data=merged_foreground, x='overlap_group', y='Mean_Essentiality_Score', 
                  color='black', alpha=0.3, size=3, jitter=0.2, order=groups)
    plt.title('Essentiality vs Data Overlap', fontsize=18, fontweight='bold', pad=20)
    plt.xlabel('Overlap Count (Number of Datasets)', fontsize=14)
    plt.ylabel('Essentiality Score (Mean Log2FC)', fontsize=14)
    format_ax(ax)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Overlap_Boxplot.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Overlap_Boxplot.pdf'), dpi=300, bbox_inches='tight')
    plt.close()

    # 6. Consistency Heatmap
    # ... [Heatmap code] ...
    top_genes_list = merged_foreground.sort_values('Mean_Essentiality_Score').head(20)['Gene'].tolist()
    heatmap_data = screen_data[screen_data['Gene'].isin(top_genes_list)]
    heatmap_matrix = heatmap_data.pivot(index='Gene', columns='CellLine', values='Score')
    heatmap_matrix['mean'] = heatmap_matrix.mean(axis=1)
    heatmap_matrix = heatmap_matrix.sort_values('mean').drop('mean', axis=1)
    
    plt.figure(figsize=(8, 10))
    sns.heatmap(heatmap_matrix, cmap='RdPu_r', annot=True, fmt='.2f', 
                linewidths=1, linecolor='white', cbar_kws={'label': 'Log2FC'})
    # 6. Consistency Heatmap
    # ... [Heatmap code] ...
    top_genes_list = merged_foreground.sort_values('Mean_Essentiality_Score').head(20)['Gene'].tolist()
    heatmap_data = screen_data[screen_data['Gene'].isin(top_genes_list)]
    heatmap_matrix = heatmap_data.pivot(index='Gene', columns='CellLine', values='Score')
    heatmap_matrix['mean'] = heatmap_matrix.mean(axis=1)
    heatmap_matrix = heatmap_matrix.sort_values('mean').drop('mean', axis=1)
    
    plt.figure(figsize=(8, 10))
    sns.heatmap(heatmap_matrix, cmap='RdPu_r', annot=True, fmt='.2f', 
                linewidths=1, linecolor='white', cbar_kws={'label': 'Log2FC'})
    plt.title('Essentiality Consistency (Top 20 Genes)', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('')
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Consistency_Heatmap.pdf'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Consistency_Heatmap.pdf'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- 7. Expression & LFC Analysis ---
    print("Generating Expression & LFC Plots...")
    
    # Load Master Ortholog Matrix
    MATRIX_FILE = "streamlit_deg_explorer/master_ortholog_matrix.csv.gz"
    print(f"Loading Master Matrix from {MATRIX_FILE}...")
    try:
        master_df = pd.read_csv(MATRIX_FILE)
        
        # Calculate Means
        # Identify columns dynamically
        tpm_cols = [c for c in master_df.columns if '_tpm' in c]
        lfc_cols = [c for c in master_df.columns if '_lfc' in c]
        
        print(f"  Found {len(tpm_cols)} TPM cols and {len(lfc_cols)} LFC cols.")
        
        master_df['Mean_TPM'] = master_df[tpm_cols].mean(axis=1)
        master_df['Mean_LFC'] = master_df[lfc_cols].mean(axis=1)
        
        # Merge with merged_foreground (Essentiality Data)
        # merged_foreground has 'Gene' (Symbol). master_df has 'Symbol' (last col from header inspection) or 'human_symbol' usually?
        # Header inspection showed 'Symbol' as the last column.
        # Let's check if 'Symbol' exists, else try 'human_symbol' or 'Human_Symbol'
        symbol_col = 'Symbol' if 'Symbol' in master_df.columns else None
        if not symbol_col:
             # Fallback
             possible = [c for c in master_df.columns if 'symbol' in c.lower()]
             if possible: symbol_col = possible[0]
        
        if symbol_col:
            print(f"  Merging on {symbol_col}...")
            # We need to merge the Expression Data onto our Core Genes
            # merged_foreground contains the Core Genes with Essentiality Scores
            expression_merged = merged_foreground.merge(master_df[[symbol_col, 'Mean_TPM', 'Mean_LFC']], 
                                                        left_on='Gene', right_on=symbol_col, how='left')
            
            # Define Groups by Species Origin
            # 1. Non-Essential: Score >= -0.5
            # 2. Mouse Specific: Score < -0.5 & Origin == 'Mouse'
            # 3. Human Specific: Score < -0.5 & Origin == 'Human'
            # 4. Shared: Score < -0.5 & Origin == 'Both' (or 'Shared')
            
            def define_group(row):
                score = row['Mean_Essentiality_Score']
                origin = str(row['Origin'])
                
                if pd.isna(score): return 'Unknown'
                if score >= -0.5: return 'Non-Essential'
                
                if origin == 'Mouse': return 'Mouse Specific'
                if origin == 'Human': return 'Human Specific'
                if origin == 'Both': return 'Shared'
                return 'Other' # Should be rare
            
            expression_merged['Group'] = expression_merged.apply(define_group, axis=1)
            
            # Filter unknowns
            plot_df = expression_merged[expression_merged['Group'] != 'Unknown'].copy()
            
            # Order
            order = ['Non-Essential', 'Mouse Specific', 'Human Specific', 'Shared']
            plot_df = plot_df[plot_df['Group'].isin(order)]
            
            # Stats (simple t-test pairwise against Non-Essential?)
            # Or just plot for now.
            print(f"  Plotting groups: {plot_df['Group'].value_counts()}")
            
            # A. TPM Plot
            plt.figure(figsize=(10, 8))
            ax = plt.gca()
            # Log transform TPM
            plot_df['Log_TPM'] = np.log2(plot_df['Mean_TPM'] + 1)
            
            sns.boxplot(data=plot_df, x='Group', y='Log_TPM', order=order, palette='viridis', width=0.5,
                        boxprops=dict(alpha=0.7))
            sns.stripplot(data=plot_df, x='Group', y='Log_TPM', order=order, color='black', alpha=0.3, size=2, jitter=0.2)
            
            plt.title('Expression vs Essentiality (By Species)', fontsize=18, fontweight='bold', pad=20)
            plt.xlabel('Essentiality Group', fontsize=14)
            plt.ylabel('Mean Expression (log2 TPM+1)', fontsize=14)
            format_ax(ax)
            plt.tight_layout()
            plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_vs_Expression_TPM.pdf'), dpi=300, bbox_inches='tight')
            plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_vs_Expression_TPM.pdf'), dpi=300, bbox_inches='tight')
            plt.close()
            plt.close()
            
            # B. LFC Plot
            plt.figure(figsize=(10, 8))
            ax = plt.gca()
            
            sns.boxplot(data=plot_df, x='Group', y='Mean_LFC', order=order, palette='viridis', width=0.5,
                         boxprops=dict(alpha=0.7))
            sns.stripplot(data=plot_df, x='Group', y='Mean_LFC', order=order, color='black', alpha=0.3, size=2, jitter=0.2)
            
            plt.title('Fold-Change vs Essentiality (By Species)', fontsize=18, fontweight='bold', pad=20)
            plt.xlabel('Essentiality Group', fontsize=14)
            plt.ylabel('Mean Log2 Fold-Change', fontsize=14)
            # Add horizontal line at 0
            plt.axhline(0, color='gray', linestyle='--')
            format_ax(ax)
            plt.tight_layout()
            plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_vs_LFC.pdf'), dpi=300, bbox_inches='tight')
            plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_vs_LFC.pdf'), dpi=300, bbox_inches='tight')
            plt.close()
            plt.close()
            
            # Save Extended Data
            print("Saving extended csv with Expression data...")
            expression_merged.to_csv(os.path.join(OUTPUT_DIR, 'Core_DEGs_Essentiality_Expression.csv'), index=False)
            
        else:
            print("  ERROR: Could not find Symbol column in Master Matrix.")
            
    except Exception as e:
        print(f"  ERROR Processing Expression Data: {e}")
        import traceback
        traceback.print_exc()

    print(f"Plots saved to {OUTPUT_DIR}")

    print(f"Plots saved to {OUTPUT_DIR}")

    # --- 8. Phase 2: DepMap Integration (Liver Context) ---
    print("\n--- Checking for DepMap Data (Phase 2) ---")
    DEPMAP_DIR = "Analysis/downstream_analysis/essentiality"
    MODEL_FILE = os.path.join(DEPMAP_DIR, "Model.csv")
    EFFECT_FILE = os.path.join(DEPMAP_DIR, "CRISPRGeneEffect.csv")
    
    if os.path.exists(MODEL_FILE) and os.path.exists(EFFECT_FILE):
        print("DepMap files found. Proceeding with Unified (Liang + DepMap) Analysis...")
        try:
            # 1. Load Model and Filter for Liver
            print("  Loading Model.csv...")
            model_df = pd.read_csv(MODEL_FILE)
            lineage_col = 'OncotreeLineage' if 'OncotreeLineage' in model_df.columns else None
            
            if lineage_col:
                # Find Liver models and mapping
                liver_models_df = model_df[model_df[lineage_col] == 'Liver']
                liver_models = liver_models_df['ModelID'].tolist()
                
                # Create ID -> Name mapping
                name_col = 'CellLineName' if 'CellLineName' in model_df.columns else 'StrippedCellLineName'
                if 'StrippedCellLineName' in model_df.columns: name_col = 'StrippedCellLineName' # Prefer stripped
                
                id_to_name = liver_models_df.set_index('ModelID')[name_col].to_dict()
                
                print(f"  Found {len(liver_models)} Liver cell lines.")
                
                if liver_models:
                    # 2. Extract Data
                    print(f"  Loading CRISPRGeneEffect.csv (cols for {len(liver_models)} liver lines)...")
                    full_effect = pd.read_csv(EFFECT_FILE, index_col=0)
                    liver_effect = full_effect[full_effect.index.isin(liver_models)]
                    
                    # RENAME INDEX to Cell Line Names
                    liver_effect.rename(index=id_to_name, inplace=True)
                    # Handle duplicate cell line names (averaging)
                    liver_effect = liver_effect.groupby(level=0).mean()
                    
                    # Clean Column Names "Symbol (ID)" -> "Symbol"
                    depmap_genes = pd.Series(liver_effect.columns)
                    depmap_symbols = depmap_genes.str.split(' ').str[0]
                    liver_effect.columns = depmap_symbols
                    
                    # Filter for duplicates
                    liver_effect = liver_effect.loc[:, ~liver_effect.columns.duplicated()]
                    
                    # Filter for Core Genes
                    # target_genes = set(merged_foreground['Gene'])  # OLD: Restricted to Liang genes
                    target_genes = foreground_symbols  # NEW: All Core Genes
                    print(f"  Target Genes for DepMap: {len(target_genes)}")
                    available_genes = list(set(liver_effect.columns) & target_genes)
                    print(f"  Found {len(available_genes)} Core Genes in DepMap.")
                    
                    liver_effect_core = liver_effect[available_genes]
                    
                    # 3. Calculate DepMap Metrics
                    depmap_mean = liver_effect_core.mean(axis=0)
                    depmap_mean.name = 'DepMap_Cas9_Score'
                    
                    # Identify Essentiality Overlap (Score < -0.5)
                    # Count number of liver lines where gene is essential
                    depmap_overlap = (liver_effect_core < -0.5).sum(axis=0)
                    depmap_overlap.name = 'DepMap_Overlap_Count'
                    
                    # ---------------------------------------------------------
                    # RE-BUILD merged_foreground to be INCLUSIVE of ALL Core Genes
                    # ---------------------------------------------------------
                    print("  Re-building merged_foreground to include ALL Core DEGs (Union of Liang & DepMap)...")
                    
                    # 1. Start with All Core Genes
                    all_genes_df = pd.DataFrame({'Gene': list(foreground_symbols)})
                    
                    # 2. Merge Metadata (Origin)
                    # meta_df has 'gene_symbol', 'Origin'
                    all_genes_df = all_genes_df.merge(meta_df[['gene_symbol', 'Origin']], 
                                                    left_on='Gene', right_on='gene_symbol', how='left')
                    
                    # 3. Merge Liang Scores (Mean Essentiality)
                    # mean_scores has 'Gene', 'Mean_Essentiality_Score'
                    all_genes_df = all_genes_df.merge(mean_scores[['Gene', 'Mean_Essentiality_Score']], 
                                                    on='Gene', how='left')
                    
                    # 4. Merge Liang Overlap Count
                    # core_df has 'gene_symbol', 'overlap_count'
                    # Need to drop duplicates because core_df might be long format? No, core_df unique on symbol?
                    # meta_df used drop_duplicates, so we can use meta_df for overlap_count too if it was kept
                    # Let's check meta_df creation: meta_df = core_df[['gene_symbol', 'overlap_count', 'analyses']]...
                    # So meta_df has overlap_count.
                    all_genes_df = all_genes_df.merge(meta_df[['gene_symbol', 'overlap_count']], 
                                                    left_on='Gene', right_on='gene_symbol', how='left', suffixes=('', '_dup'))
                    # Clean up
                    if 'overlap_count_dup' in all_genes_df.columns:
                        all_genes_df['overlap_count'] = all_genes_df['overlap_count'].fillna(all_genes_df['overlap_count_dup'])
                        all_genes_df.drop(columns=['overlap_count_dup', 'gene_symbol', 'gene_symbol_dup'], errors='ignore', inplace=True)
                    else:
                        all_genes_df.drop(columns=['gene_symbol'], errors='ignore', inplace=True)

                    # 5. Merge DepMap Scores (DepMap_Cas9_Score)
                    all_genes_df = all_genes_df.merge(depmap_mean, left_on='Gene', right_index=True, how='left')
                    all_genes_df = all_genes_df.merge(depmap_overlap, left_on='Gene', right_index=True, how='left')
                    
                    # Update the global variable
                    merged_foreground = all_genes_df
                    
                    # ---------------------------------------------------------
                    # 4. UNIFIED PLOTTING
                    # ---------------------------------------------------------
                    print("  Generating Unified Plots...")
                    
                    # A. Unified Density (Liang vs DepMap)
                    # Prepare long format
                    liang_scores = merged_foreground[['Mean_Essentiality_Score']].dropna()
                    liang_scores['Method'] = 'Liang (Cas13 RNA-KD)'
                    liang_scores.rename(columns={'Mean_Essentiality_Score': 'Score'}, inplace=True)
                    
                    depmap_scores = merged_foreground[['DepMap_Cas9_Score']].dropna()
                    depmap_scores['Method'] = 'DepMap (Cas9 DNA-KO)'
                    depmap_scores.rename(columns={'DepMap_Cas9_Score': 'Score'}, inplace=True)
                    
                    unified_density = pd.concat([liang_scores, depmap_scores])
                    
                    plt.figure(figsize=(10, 6))
                    ax = plt.gca()
                    sns.kdeplot(data=unified_density, x='Score', hue='Method', fill=True, palette=['#d358c7', 'orange'], alpha=0.4)
                    plt.title('Distribution of Essentiality: Liang (RNA) vs DepMap (DNA)', fontsize=16)
                    plt.xlabel('Essentiality Score', fontsize=14)
                    plt.axvline(-0.5, linestyle='--', color='gray', label='Essential Threshold')
                    format_ax(ax)
                    plt.tight_layout()
                    plt.savefig(os.path.join(OUTPUT_DIR, 'Unified_Density_Liang_vs_DepMap.pdf'), dpi=300, bbox_inches='tight')
                    plt.savefig(os.path.join(OUTPUT_DIR, 'Unified_Density_Liang_vs_DepMap.pdf'), dpi=300, bbox_inches='tight')
                    plt.close()

                    # B. Unified Scatter (Comparison) with Marginal Strip for Missing Data
                    print("  Generating Unified Scatter with Marginals...")
                    from matplotlib.gridspec import GridSpec
                    
                    # Data Preparation
                    df = merged_foreground.copy()
                    ESSENTIAL_THRESH = -0.5
                    
                    # Define Categories for Coloring
                    def get_scatter_category(row):
                         l = row['Mean_Essentiality_Score']
                         d = row['DepMap_Cas9_Score']
                         
                         # Check strict NaNs
                         l_missing = pd.isna(l)
                         d_missing = pd.isna(d)
                         
                         if l_missing and d_missing: return 'Missing Data'
                         if l_missing: return 'DepMap Only (Not Screened in Liang)' # Plotted on Marginal
                         if d_missing: return 'Liang Only (Not in DepMap)'         # Not plotted (rare)
                         
                         # Both Present - Determine Essentiality Intersection
                         l_ess = l < ESSENTIAL_THRESH
                         d_ess = d < ESSENTIAL_THRESH
                         
                         if l_ess and d_ess: return 'Essential (Both)'
                         if l_ess: return 'Liang Specific Essential'
                         if d_ess: return 'DepMap Specific Essential'
                         return 'Non-Essential (Both)'

                    df['Category'] = df.apply(get_scatter_category, axis=1)
                    
                    # Separate Data
                    scatter_df = df[~df['Category'].str.contains('Only') & (df['Category'] != 'Missing Data')].copy()
                    depmap_only_df = df[df['Category'] == 'DepMap Only (Not Screened in Liang)'].copy()
                    
                    # Calculate Counts for Legend
                    counts = scatter_df['Category'].value_counts()
                    
                    # Sort scatter_df to put 'Essential (Both)' on top (z-order)
                    cat_order = ['Non-Essential (Both)', 'Liang Specific Essential', 'DepMap Specific Essential', 'Essential (Both)']
                    scatter_df['Category'] = pd.Categorical(scatter_df['Category'], categories=cat_order, ordered=True)
                    scatter_df = scatter_df.sort_values('Category')
                    
                    # Rename Categories to Include Counts (and simplify labels)
                    # Mapping: Internal Name -> Display Name
                    name_map = {
                        'Essential (Both)': 'Both',
                        'Liang Specific Essential': 'Liang Only',
                        'DepMap Specific Essential': 'DepMap Only',
                        'Non-Essential (Both)': 'Non-Essential'
                    }
                    
                    cat_map = {}
                    for c in cat_order:
                        display_name = name_map.get(c, c)
                        count = counts.get(c, 0)
                        cat_map[c] = f"{display_name} (n={count})"
                        
                    scatter_df['Legend_Label'] = scatter_df['Category'].map(cat_map)
                    
                    # Colors (Map new labels to colors)
                    base_palette = {
                        'Both': '#e14b9d',          
                        'Liang Only': '#FBC02D',  
                        'DepMap Only': '#FF7043', 
                        'Non-Essential': '#E0E0E0'       
                    }
                    # Map loop to use correct keys
                    palette = {}
                    for original_cat, new_label in cat_map.items():
                         # Find which base color key matches the original category logic
                         simple_name = name_map.get(original_cat)
                         if simple_name in base_palette:
                             palette[new_label] = base_palette[simple_name]
                    
                    # Plot Setup
                    fig = plt.figure(figsize=(12, 8)) # Slightly wider
                    gs = GridSpec(1, 2, width_ratios=[4, 1.2], wspace=0.08)
                    
                    # 1. Main Scatter
                    ax_main = fig.add_subplot(gs[0])
                    sns.scatterplot(data=scatter_df, x='Mean_Essentiality_Score', y='DepMap_Cas9_Score', 
                                    hue='Legend_Label', palette=palette, alpha=0.8, s=30, ax=ax_main, edgecolor='w', linewidth=0.3)
                    
                    # Add Quadrant Counts - Explicitly for User Clarity
                    q_counts = df['Category'].value_counts()
                    # Mapping our categories to quadrants
                    # Essential (Both) -> Bottom Left
                    # Liang Specific Essential -> Bottom Right (Liang < -0.5, DepMap > -0.5) WAIT
                    # Liang Ess: Liang < -0.5. If DepMap > -0.5, then it is Liang Specific.
                    # Coordinates: x < -0.5, y > -0.5. Top Left? No.
                    # x (Liang) < -0.5 (Left). y (DepMap) > -0.5 (Top). So Top Left is Liang Specific.
                    # DepMap Specific: x > -0.5 (Right). y < -0.5 (Bottom). So Bottom Right is DepMap Specific.
                    # Non-Essential: Top Right.
                    # Essential Both: Bottom Left.
                    
                    # Let's count properly based on quadrant logic
                    valid_df = df.dropna(subset=['Mean_Essentiality_Score', 'DepMap_Cas9_Score'])
                    n_both = ((valid_df['Mean_Essentiality_Score'] < -0.5) & (valid_df['DepMap_Cas9_Score'] < -0.5)).sum()
                    n_liang_only = ((valid_df['Mean_Essentiality_Score'] < -0.5) & (valid_df['DepMap_Cas9_Score'] >= -0.5)).sum()
                    n_depmap_only = ((valid_df['Mean_Essentiality_Score'] >= -0.5) & (valid_df['DepMap_Cas9_Score'] < -0.5)).sum()
                    n_non = ((valid_df['Mean_Essentiality_Score'] >= -0.5) & (valid_df['DepMap_Cas9_Score'] >= -0.5)).sum()
                    

                    
                    # Guidelines
                    ax_main.axhline(ESSENTIAL_THRESH, linestyle='--', color='gray', alpha=0.5)
                    ax_main.axvline(ESSENTIAL_THRESH, linestyle='--', color='gray', alpha=0.5)
                    
                    # Labels
                    ax_main.set_title('Cross-Modality Essentiality Comparison', fontsize=18, fontweight='bold')
                    ax_main.set_xlabel('Liang Cas13 Score (Mean Log2FC)', fontsize=14)
                    ax_main.set_ylabel('DepMap Cas9 Score (Mean)', fontsize=14)
                    
                    # Legend (Enlarged and Clearer)
                    ax_main.legend(fontsize=12, loc='upper left', frameon=True, title='Shared Genes Status', 
                                   title_fontsize=13, framealpha=0.95, borderpad=1)
                    format_ax(ax_main)
                    
                    # 2. Right Marginal: DepMap-only Genes (Missing in Liang)
                    ax_right = fig.add_subplot(gs[1], sharey=ax_main)
                    
                    if len(depmap_only_df) > 0:
                        # Categorize DepMap Only genes by Essentiality
                        depmap_only_df['Status'] = depmap_only_df['DepMap_Cas9_Score'].apply(lambda x: 'Essential' if x < ESSENTIAL_THRESH else 'Non-Essential')
                        
                        # Calculate counts for strip plot
                        d_counts = depmap_only_df['Status'].value_counts()
                        n_ess = d_counts.get('Essential', 0)
                        n_non = d_counts.get('Non-Essential', 0)
                        
                        # Split Data for Layering
                        depmap_only_ess = depmap_only_df[depmap_only_df['Status'] == 'Essential']
                        depmap_only_non = depmap_only_df[depmap_only_df['Status'] == 'Non-Essential']
                        
                        # 1. Plot Non-Essential (Background) - Gray
                        sns.stripplot(data=depmap_only_non, y='DepMap_Cas9_Score', color='#BDBDBD',
                                      alpha=0.3, s=3, jitter=0.4, ax=ax_right, zorder=1)
                        
                        # 2. Plot Essential (Foreground) - Orange
                        # Higher alpha, zorder=2 ensures they are ON TOP
                        sns.stripplot(data=depmap_only_ess, y='DepMap_Cas9_Score', color='#FF7043',
                                      alpha=0.9, s=4.5, jitter=0.4, ax=ax_right, zorder=2)
                                      
                        # Annotate
                        n_miss = len(depmap_only_df)
                        ax_right.set_xlabel(f'Not Screened\nin Liang\n(n={n_miss})', fontsize=12, fontweight='bold')
                        ax_right.axhline(ESSENTIAL_THRESH, linestyle='--', color='gray', alpha=0.5)
                        
                        # Add explicit text labels for counts
                        # Essential - Moved inside/closer, Black color
                        # Using transAxes to place it at the bottom center of the panel (where essential genes are)
                        ax_right.text(0.5, 0.05, f"n={n_ess}", ha='center', va='bottom', 
                                      fontsize=12, color='black', fontweight='bold', transform=ax_right.transAxes,
                                      bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=2))
                        
                        # Removed Non-Essential label as requested

                    else:
                        ax_right.text(0.5, 0.5, "No Missing\nLiang Data", ha='center')
                        
                    ax_right.set_ylabel('')
                    ax_right.tick_params(labelleft=False) # Hide Y ticks
                    sns.despine(ax=ax_right, left=True, bottom=False)
                    ax_right.set_facecolor('#fdfdfd') # Very light gray options
                    
                    plt.suptitle('Unified Essentiality Landscape', fontsize=20, fontweight='bold', y=0.98)
                    plt.savefig(os.path.join(OUTPUT_DIR, 'Unified_Scatter_Comparison.pdf'), dpi=300, bbox_inches='tight')
                    plt.savefig(os.path.join(OUTPUT_DIR, 'Unified_Scatter_Comparison.pdf'), dpi=300, bbox_inches='tight')
                    plt.close()

                    # C. Unified Heatmap
                    # Pick Top 20 Essential Genes (based on Liang score, to show consistency across both)
                    top_genes = merged_foreground.sort_values('Mean_Essentiality_Score').drop_duplicates('Gene').head(20)['Gene'].tolist()
                    print(f"DEBUG: Top Genes: {len(top_genes)} (Unique: {len(set(top_genes))})")
                    
                    # Get Data: Liang Lines (need to re-merge or extract if available)
                    # Note: We aggregated Liang data earlier. If we want raw lines, we need 'data_matrix'.
                    # Assuming 'data_matrix' from Step 4 is still available.
                    
                    # Prepare Heatmap Data
                    # Liang Matrix (Top 20)
                    # Reconstruct from screen_data which is available globally
                    liang_data_subset = screen_data[screen_data['Gene'].isin(top_genes)]
                    liang_data_subset = liang_data_subset.drop_duplicates(subset=['Gene', 'CellLine'])
                    liang_heatmap = liang_data_subset.pivot(index='Gene', columns='CellLine', values='Score')
                    print(f"DEBUG: Liang Heatmap Shape: {liang_heatmap.shape}")
                    print(f"DEBUG: Liang Index Unique? {liang_heatmap.index.is_unique}")
                    
                    # DepMap Matrix (Top 20, limited to 5 random liver lines to fit?) or Average?
                    # User said "show all cell line data together".
                    # If DepMap has 25 lines + 5 Liang = 30 cols. Doable.
                    depmap_heatmap = liver_effect_core[top_genes].T # Transpose: Index=Gene, Col=Model
                    print(f"DEBUG: DepMap Heatmap Shape: {depmap_heatmap.shape}")
                    print(f"DEBUG: DepMap Heatmap Index Unique? {depmap_heatmap.index.is_unique}")
                    
                    # Combine
                    # Rename Liang cols to "Liang_..."
                    liang_renamed = liang_heatmap.add_prefix('Liang_')
                    depmap_renamed = depmap_heatmap.add_prefix('DepMap_')
                    
                    # Combine using Merge (safer)
                    df1 = liang_renamed.reset_index()
                    df2 = depmap_renamed.reset_index()
                    # Ensure 'Gene' col exists in df2 (it was index)
                    if 'index' in df2.columns: df2.rename(columns={'index': 'Gene'}, inplace=True)
                    
                    combined_heatmap = pd.merge(df1, df2, on='Gene', how='outer').set_index('Gene')
                    
                    if not combined_heatmap.index.is_unique:
                        print("Reducing duplicates in Combined Heatmap...")
                        combined_heatmap = combined_heatmap.groupby(level=0).mean()
                    
                    # Save Data for inspection
                    combined_heatmap.to_csv(os.path.join(OUTPUT_DIR, 'Unified_Heatmap_Data.csv'))
                    
                    try:
                        plt.figure(figsize=(12, 10))
                        sns.heatmap(combined_heatmap, cmap='RdBu', center=0, cbar_kws={'label': 'Essentiality Score'})
                        plt.title('Consistency of Top Essential Genes (Liang vs DepMap)', fontsize=16)
                        plt.tight_layout()
                        plt.savefig(os.path.join(OUTPUT_DIR, 'Unified_Heatmap.pdf'), dpi=300, bbox_inches='tight')
                        plt.savefig(os.path.join(OUTPUT_DIR, 'Unified_Heatmap.pdf'), dpi=300, bbox_inches='tight')
                        plt.close()
                    except Exception as e:
                        print(f"ERROR plotting Unified Heatmap: {e}")
                    
                    # D. Unified Overlap plot (Percent Essential)
                    try:
                        # Problem: Liang has 5 lines, DepMap has N. "Count" is not comparable. "Percentage" is.
                        merged_foreground['Liang_Pct_Essential'] = merged_foreground['overlap_count'] / 5 * 100
                        merged_foreground['DepMap_Pct_Essential'] = merged_foreground['DepMap_Overlap_Count'] / len(liver_models) * 100
                        
                        # Melt for boxplot
                        mp1 = merged_foreground[['Liang_Pct_Essential', 'Mean_Essentiality_Score']].copy()
                        mp1.columns = ['Pct_Essential', 'Score']
                        mp1['Dataset'] = 'Liang (Cas13)'
                        
                        mp2 = merged_foreground[['DepMap_Pct_Essential', 'DepMap_Cas9_Score']].copy()
                        mp2.columns = ['Pct_Essential', 'Score']
                        mp2['Dataset'] = 'DepMap (Cas9)'
                        
                        unified_overlap = pd.concat([mp1, mp2]).reset_index(drop=True)
                        
                        # Plot: Dist of Scores for "Shared Essential" (>80% lines) vs "Specific" (<20%)
                        # Define categories
                        def get_status(pct):
                            if pct > 80: return 'Shared (>80%)'
                            if pct < 20: return 'Rare (<20%)'
                            return 'Partial'
                            
                        unified_overlap['Consistency'] = unified_overlap['Pct_Essential'].apply(get_status)
                        
                        plt.figure(figsize=(10,6))
                        ax = plt.gca()
                        sns.boxplot(data=unified_overlap, x='Consistency', y='Score', hue='Dataset', palette=['#d358c7', 'orange'], order=['Rare (<20%)', 'Partial', 'Shared (>80%)'])
                        plt.title('Essentiality Score by Consistency (Pct Lines Essential)', fontsize=16)
                        plt.axhline(-0.5, linestyle='--', color='gray')
                        format_ax(ax)
                        plt.tight_layout()
                        plt.savefig(os.path.join(OUTPUT_DIR, 'Unified_Consistency_Boxplot.pdf'), dpi=300, bbox_inches='tight')
                        plt.savefig(os.path.join(OUTPUT_DIR, 'Unified_Consistency_Boxplot.pdf'), dpi=300, bbox_inches='tight')
                        plt.close()
                    except Exception as e:
                        print(f"ERROR plotting Unified Boxplot: {e}")

                    print("  Saved Unified Plots.")
                    merged_foreground.to_csv(os.path.join(OUTPUT_DIR, 'Unified_Essentiality_Combined.csv'), index=False)
                    
                    # ---------------------------------------------------------
                    # 5. REGENERATED COMPARISON PLOTS (Side-by-Side Liang vs DepMap)
                    # ---------------------------------------------------------
                    print("  Generating Side-by-Side Comparison Plots...")
                    
                    # Prepare DepMap data for comparison plots
                    depmap_scores_df = merged_foreground[['Gene', 'DepMap_Cas9_Score', 'DepMap_Overlap_Count']].dropna()
                    
                    # A. Density Comparison (Overlaid with Both Backgrounds + Core)
                    plt.figure(figsize=(10, 7))
                    ax = plt.gca()
                    
                    # 1. Liang Background (All genes in Screen not in Core)
                    # background_in_screen is defined earlier
                    sns.kdeplot(data=background_in_screen, x='Mean_Essentiality_Score', 
                                label=f'Liang Background (n={len(background_in_screen)})', 
                                color='#bdbdbd', fill=True, alpha=0.2, linewidth=0)
                    
                    # 2. DepMap Background (All Liver genes not in Core)
                    # We need to get all liver genes, exclude Core
                    all_liver_genes = set(liver_effect.columns)
                    depmap_bg_genes = list(all_liver_genes - target_genes)
                    if len(depmap_bg_genes) > 0:
                         depmap_bg_metrics = liver_effect[depmap_bg_genes].mean(axis=0).to_frame(name='Score')
                         sns.kdeplot(data=depmap_bg_metrics, x='Score',
                                     label=f'DepMap Background (n={len(depmap_bg_genes)})',
                                     color='#ffcc80', fill=True, alpha=0.2, linewidth=0) # Light Orange

                    # 3. Liang Core
                    sns.kdeplot(data=foreground_in_screen, x='Mean_Essentiality_Score', 
                                label=f'Liang Core DEGs (n={len(foreground_in_screen)})', 
                                color='#e14b9d', fill=True, alpha=0.5, linewidth=2) # Pink
                    
                    # 4. DepMap Core
                    sns.kdeplot(data=depmap_scores_df, x='DepMap_Cas9_Score', 
                                label=f'DepMap Core DEGs (n={len(depmap_scores_df)})', 
                                color='orange', fill=True, alpha=0.5, linewidth=2) # Orange
                    
                    plt.axvline(-0.5, linestyle='--', color='gray', label='Essential Threshold')
                    plt.title('Essentiality Density: Liang vs DepMap', fontsize=18, fontweight='bold')
                    plt.xlabel('Essentiality Score', fontsize=14)
                    plt.ylabel('Density', fontsize=14)
                    plt.legend(frameon=False, fontsize=10, loc='upper left')
                    format_ax(ax)
                    plt.tight_layout()
                    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Density.pdf'), dpi=300, bbox_inches='tight')
                    plt.close()

                    
                    # C. Raincloud/Boxplot Comparison (Side-by-Side)
                    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)
                    
                    # Liang Raincloud
                    ax1 = axes[0]
                    df_liang = pd.DataFrame({'Score': foreground_in_screen['Mean_Essentiality_Score'], 'Group': 'Core DEGs'})
                    sns.boxplot(data=df_liang, x='Group', y='Score', color=COLOR_CORE, width=0.4, ax=ax1)
                    sns.stripplot(data=df_liang, x='Group', y='Score', color='black', alpha=0.3, size=2, ax=ax1)
                    ax1.axhline(-0.5, linestyle='--', color='gray', alpha=0.7)
                    ax1.set_title('Liang (Cas13)', fontsize=14, fontweight='bold')
                    ax1.set_xlabel('')
                    ax1.set_ylabel('Essentiality Score', fontsize=12)
                    format_ax(ax1)
                    
                    # DepMap Raincloud
                    ax2 = axes[1]
                    df_depmap = pd.DataFrame({'Score': depmap_scores_df['DepMap_Cas9_Score'], 'Group': 'Core DEGs'})
                    sns.boxplot(data=df_depmap, x='Group', y='Score', color='orange', width=0.4, ax=ax2)
                    sns.stripplot(data=df_depmap, x='Group', y='Score', color='black', alpha=0.3, size=2, ax=ax2)
                    ax2.axhline(-0.5, linestyle='--', color='gray', alpha=0.7)
                    ax2.set_title('DepMap (Cas9)', fontsize=14, fontweight='bold')
                    ax2.set_xlabel('')
                    format_ax(ax2)
                    
                    plt.suptitle('Essentiality Distribution: Cas13 vs Cas9', fontsize=16, fontweight='bold', y=1.02)
                    plt.tight_layout()
                    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Raincloud.pdf'), dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    # D. Overlap Boxplot Comparison (Side-by-Side)
                    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
                    
                    # Liang Overlap (1-5 cell lines)
                    ax1 = axes[0]
                    merged_foreground['Liang_overlap_group'] = merged_foreground['overlap_count'].apply(
                        lambda x: str(int(x)) if pd.notnull(x) else 'Unknown')
                    liang_groups = sorted([g for g in merged_foreground['Liang_overlap_group'].unique() if g != 'Unknown'])
                    sns.boxplot(data=merged_foreground[merged_foreground['Liang_overlap_group'] != 'Unknown'], 
                                x='Liang_overlap_group', y='Mean_Essentiality_Score', 
                                palette='magma_r', order=liang_groups, width=0.5, ax=ax1)
                    ax1.axhline(-0.5, linestyle='--', color='gray', alpha=0.7)
                    ax1.set_title('Liang (Cas13) - 5 Cell Lines', fontsize=14, fontweight='bold')
                    ax1.set_xlabel('# Cell Lines Essential', fontsize=12)
                    ax1.set_ylabel('Essentiality Score', fontsize=12)
                    format_ax(ax1)
                    
                    # DepMap Overlap (Binned: 0-10%, 10-50%, 50-90%, 90-100%)
                    ax2 = axes[1]
                    n_liver = len(liver_models)
                    merged_foreground['DepMap_overlap_pct'] = merged_foreground['DepMap_Overlap_Count'] / n_liver * 100
                    def bin_depmap_overlap(pct):
                        if pd.isna(pct): return 'N/A'
                        if pct < 10: return '<10%'
                        if pct < 50: return '10-50%'
                        if pct < 90: return '50-90%'
                        return '>90%'
                    merged_foreground['DepMap_overlap_group'] = merged_foreground['DepMap_overlap_pct'].apply(bin_depmap_overlap)
                    depmap_order = ['<10%', '10-50%', '50-90%', '>90%']
                    depmap_plot_df = merged_foreground[merged_foreground['DepMap_overlap_group'] != 'N/A']
                    sns.boxplot(data=depmap_plot_df, x='DepMap_overlap_group', y='DepMap_Cas9_Score', 
                                palette='YlOrRd', order=depmap_order, width=0.5, ax=ax2)
                    ax2.axhline(-0.5, linestyle='--', color='gray', alpha=0.7)
                    ax2.set_title(f'DepMap (Cas9) - {n_liver} Liver Lines', fontsize=14, fontweight='bold')
                    ax2.set_xlabel('% Cell Lines Essential', fontsize=12)
                    format_ax(ax2)
                    
                    plt.suptitle('Essentiality vs Consistency: Liang vs DepMap', fontsize=16, fontweight='bold', y=1.02)
                    plt.tight_layout()
                    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Overlap_Boxplot.pdf'), dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    # E. Species-Split Violin (Side-by-Side)
                    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
                    
                    ax1 = axes[0]
                    sns.violinplot(data=merged_foreground, x='Origin', y='Mean_Essentiality_Score', 
                                   hue='Origin', palette='viridis', cut=0, inner='quartile', ax=ax1, legend=False)
                    ax1.axhline(-0.5, linestyle='--', color='gray', alpha=0.7)
                    ax1.set_title('Liang (Cas13)', fontsize=14, fontweight='bold')
                    ax1.set_xlabel('Gene Origin', fontsize=12)
                    ax1.set_ylabel('Essentiality Score', fontsize=12)
                    format_ax(ax1)
                    
                    ax2 = axes[1]
                    depmap_violin_df = merged_foreground.dropna(subset=['DepMap_Cas9_Score'])
                    sns.violinplot(data=depmap_violin_df, x='Origin', y='DepMap_Cas9_Score', 
                                   hue='Origin', palette='viridis', cut=0, inner='quartile', ax=ax2, legend=False)
                    ax2.axhline(-0.5, linestyle='--', color='gray', alpha=0.7)
                    ax2.set_title('DepMap (Cas9)', fontsize=14, fontweight='bold')
                    ax2.set_xlabel('Gene Origin', fontsize=12)
                    format_ax(ax2)
                    
                    plt.suptitle('Essentiality by Species Origin: Liang vs DepMap', fontsize=16, fontweight='bold', y=1.02)
                    plt.tight_layout()
                    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Species_Split_Violin.pdf'), dpi=300, bbox_inches='tight')
                    plt.close()
                    
                    # F. Consistency Heatmap (Combined: Liang + DepMap columns)
                    # This is already done as Unified_Heatmap.png, copy it
                    # Or generate a new Essentiality_Consistency_Heatmap with same logic
                    # Just rename the file for consistency
                    import shutil
                    src_hm = os.path.join(OUTPUT_DIR, 'Unified_Heatmap.pdf')
                    dst_hm = os.path.join(OUTPUT_DIR, 'Essentiality_Consistency_Heatmap.pdf')
                    if os.path.exists(src_hm):
                        shutil.copy(src_hm, dst_hm)
                    
                    # G. Overview Stacked Bar (Add DepMap Essentiality Status)
                    # Recalculate with DepMap status
                    all_core_meta = meta_df[['gene_symbol', 'Origin']].copy()
                    all_core_meta['Liang_InScreen'] = all_core_meta['gene_symbol'].isin(foreground_in_screen['Gene'])
                    
                    # Merge Liang Scores (from merged_foreground)
                    all_core_meta = all_core_meta.merge(merged_foreground[['Gene', 'Mean_Essentiality_Score']], 
                                                        left_on='gene_symbol', right_on='Gene', how='left')
                    
                    # Merge DepMap Scores (directly from source to avoid Liang-restriction)
                    if 'depmap_mean' in locals() and depmap_mean is not None:
                        all_core_meta = all_core_meta.merge(depmap_mean, left_on='gene_symbol', right_index=True, how='left')
                    elif 'DepMap_Cas9_Score' in merged_foreground.columns:
                        # Fallback if depmap_mean var not found but column exists (should not happen with new logic)
                        all_core_meta = all_core_meta.merge(merged_foreground[['Gene', 'DepMap_Cas9_Score']], 
                                                            left_on='gene_symbol', right_on='Gene', how='left')
                    
                    ESSENTIAL_THRESH = -0.5
                    all_core_meta['Liang_Status'] = 'Not screened'
                    all_core_meta.loc[all_core_meta['Liang_InScreen'], 'Liang_Status'] = 'Non-Essential'
                    all_core_meta.loc[(all_core_meta['Liang_InScreen']) & (all_core_meta['Mean_Essentiality_Score'] < ESSENTIAL_THRESH), 'Liang_Status'] = 'Essential'
                    
                    all_core_meta['DepMap_Status'] = 'Not screened'
                    all_core_meta.loc[all_core_meta['DepMap_Cas9_Score'].notna(), 'DepMap_Status'] = 'Non-Essential'
                    all_core_meta.loc[all_core_meta['DepMap_Cas9_Score'] < ESSENTIAL_THRESH, 'DepMap_Status'] = 'Essential'
                    
                    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
                    
                    # Colors
                    colors_liang = ['#e14b9d', '#9c9c9c', '#e6e6e5']  # Pink, Gray, LightGray
                    colors_depmap = ['orange', '#9c9c9c', '#e6e6e5']  # Orange, Gray, LightGray
                    cols = ['Essential', 'Non-Essential', 'Not screened']
                    
                    # 1. Liang Plot - Simple Stacked Bar
                    ax1 = axes[0]
                    pivot_liang = all_core_meta.groupby(['Origin', 'Liang_Status']).size().unstack(fill_value=0)
                    for c in cols:
                        if c not in pivot_liang.columns: pivot_liang[c] = 0
                    pivot_liang = pivot_liang[cols]
                    pivot_liang.plot(kind='bar', stacked=True, ax=ax1, color=colors_liang, width=0.7)
                    ax1.set_title('Liang (Cas13 RNA-KD)', fontsize=18, fontweight='bold')
                    ax1.set_ylabel('Number of Genes', fontsize=16)
                    ax1.set_xlabel('Gene Origin', fontsize=16)
                    ax1.set_xticklabels(pivot_liang.index, rotation=0, fontsize=14)
                    ax1.legend(title='Status', frameon=False, fontsize=12, title_fontsize=14)
                    format_ax(ax1)
                    
                    # 2. DepMap Plot - Simple Stacked Bar
                    ax2 = axes[1]
                    pivot_depmap = all_core_meta.groupby(['Origin', 'DepMap_Status']).size().unstack(fill_value=0)
                    for c in cols:
                        if c not in pivot_depmap.columns: pivot_depmap[c] = 0
                    pivot_depmap = pivot_depmap[cols]
                    pivot_depmap.plot(kind='bar', stacked=True, ax=ax2, color=colors_depmap, width=0.7)
                    ax2.set_title('DepMap (Cas9 DNA-KO)', fontsize=18, fontweight='bold')
                    ax2.set_ylabel('Number of Genes', fontsize=16)
                    ax2.set_xlabel('Gene Origin', fontsize=16)
                    ax2.set_xticklabels(pivot_depmap.index, rotation=0, fontsize=14)
                    ax2.legend(title='Status', frameon=False, fontsize=12, title_fontsize=14)
                    format_ax(ax2)
                    
                    plt.suptitle('Essentiality Overview by Origin: Liang vs DepMap', fontsize=20, fontweight='bold', y=1.02)
                    plt.tight_layout()
                    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Overview_StackedBar.pdf'), dpi=300, bbox_inches='tight')
                    plt.savefig(os.path.join(OUTPUT_DIR, 'Essentiality_Overview_StackedBar.pdf'), dpi=300, bbox_inches='tight')
                    plt.close()


                    
                    print("  Saved Comparison Plots.")
                    
            else:
                 print("  WARNING: Could not find Lineage column in Model.csv")
                 
        except Exception as e:
            print(f"  ERROR in Unified Analysis: {e}")
            import traceback
            traceback.print_exc()
            
    else:
        print(f"DepMap files not found in {DEPMAP_DIR}. Skipping Phase 2.")
        print("To run Phase 2: Create 'RNA-seq/depmap_data' and upload 'Model.csv' + 'CRISPR_gene_effect.csv'.")

    # Save summary stats (Base)
    print("Saving summary CSV...")
    merged_foreground.to_csv(os.path.join(OUTPUT_DIR, 'Core_DEGs_Essentiality.csv'), index=False)

print("Analysis Complete.")
