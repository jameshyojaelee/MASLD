#!/usr/bin/env python3
"""
Compare published/per-study DEGs from each MASLD dataset against the
MATCHING mega-analysis and meta-analysis contrast. PDF-only output.

Contrast matching (verified against config/human_datasets.yaml):
  Disease vs Control  → dream_results.csv / meta_analysis_results.csv
    GSE126848:  (NAFL+NASH)/2 - (Control+Control_Obese)/2
    GSE135251:  (NAFL+NASH+NASH_Fibrosis)/3 - Control
    GSE130970:  (F1+F2F3+F4)/3 - Control
    GSE213621:  (F0F1+F2+F3F4)/3 - Control
    GSE162694:  (NASH_F0+...+NASH_F4)/5 - Control

  PRJNA512027 (Gerhard 2018) is excluded from the cohort presentation
  (L0/S0 library-prep batch is perfectly confounded with diagnosis: all
  34 controls L0, all 102 NASH S0). The pipeline still loads it in the
  fibrosis-vs-healthy contrast, but it is not shown here.

  NASH vs NAFL → nafl_vs_nash_dream.csv / nafl_vs_nash_meta.csv
    GSE167523:  NASH - NAFL (no healthy controls)

  Adv vs Early Fib → adv_vs_early_fibrosis_dream.csv / ..._meta.csv
    GSE174478:  (F3+F4)/2 - (F0+F1)/2  [script 15b per-study, F2 dropped]
    GSE193066:  (F3+F4)/2 - (F0+F1)/2
    GSE240729:  (F3+F4)/2 - (F0+F1)/2
"""

import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Paths
# ============================================================
BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT  = f"{BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results"
PUB_DEGS = f"{BASE}/data/published_degs"
PER_STUDY = f"{INT}/per_study"
DISEASE_SIG = f"{INT}/disease_signatures"
OUTDIR = f"{BASE}/figures"

PADJ_THRESH = 0.1
PERSTUDY_PADJ = 0.05

# ============================================================
# Load all three mega-analysis + meta-analysis references
# ============================================================
def load_mega_meta(mega_path, meta_path, mega_padj_col, mega_lfc_col,
                   meta_padj_col, meta_lfc_col):
    """Load a mega+meta pair, returning DEG sets, LFC dicts, and symbol map."""
    mega = pd.read_csv(mega_path)
    mega['gene_base'] = mega['gene'].str.replace(r'\.\d+$', '', regex=True)
    mega_degs = set(mega.loc[mega[mega_padj_col] < PADJ_THRESH, 'gene_base'])
    mega_lfc = dict(zip(mega['gene_base'], mega[mega_lfc_col]))

    sym_map = {}
    if 'symbol' in mega.columns:
        for _, row in mega.iterrows():
            if pd.notna(row['symbol']) and row['symbol'] != '':
                sym_map[row['symbol']] = row['gene_base']
    sym_lfc = {k: mega_lfc[v] for k, v in sym_map.items() if v in mega_lfc}

    meta = pd.read_csv(meta_path)
    meta['gene_base'] = meta['gene'].str.replace(r'\.\d+$', '', regex=True)
    meta_degs = set(meta.loc[meta[meta_padj_col] < PADJ_THRESH, 'gene_base'])
    meta_lfc = dict(zip(meta['gene_base'], meta[meta_lfc_col]))

    return {
        'mega_degs': mega_degs, 'mega_lfc': mega_lfc,
        'meta_degs': meta_degs, 'meta_lfc': meta_lfc,
        'sym_map': sym_map, 'sym_lfc': sym_lfc,
        'n_mega': len(mega_degs), 'n_meta': len(meta_degs),
    }

print("Loading reference contrasts...")

ref_dvc = load_mega_meta(
    f"{INT}/integration/dream_results.csv",
    f"{INT}/integration/meta_analysis_results.csv",
    'padj', 'logFC', 'meta_padj', 'meta_logFC')
print(f"  Disease vs Control: mega={ref_dvc['n_mega']:,}, meta={ref_dvc['n_meta']:,}")

ref_nn = load_mega_meta(
    f"{DISEASE_SIG}/nafl_vs_nash_dream.csv",
    f"{DISEASE_SIG}/nafl_vs_nash_meta.csv",
    'padj', 'logFC', 'meta_padj', 'meta_logFC')  # this file's padj col = 'padj' (adv_vs_early below genuinely uses 'adj.P.Val')
print(f"  NAFL vs NASH: mega={ref_nn['n_mega']:,}, meta={ref_nn['n_meta']:,}")

ref_fib = load_mega_meta(
    f"{DISEASE_SIG}/adv_vs_early_fibrosis_dream.csv",
    f"{DISEASE_SIG}/adv_vs_early_fibrosis_meta.csv",
    'adj.P.Val', 'logFC', 'meta_padj', 'meta_logFC')
print(f"  Adv vs Early Fib: mega={ref_fib['n_mega']:,}, meta={ref_fib['n_meta']:,}")

# ============================================================
# Helper: resolve gene symbols to Ensembl base IDs
# ============================================================
def symbols_to_ensembl(symbols, sym_map):
    result = set()
    upper_map = {str(s).upper(): e for s, e in sym_map.items()}
    for s in symbols:
        s_str = str(s).strip()
        if s_str in sym_map:
            result.add(sym_map[s_str])
        elif s_str.upper() in upper_map:
            result.add(upper_map[s_str.upper()])
    return result

# ============================================================
# Load published DEGs per dataset + assign correct reference
# ============================================================
print("\nLoading published DEG lists...")
datasets = {}

# --- Disease vs Control datasets (6) ---

print("  GSE130970 (Hoang)...")
gse130970_df = pd.read_excel(f"{PUB_DEGS}/GSE130970/MOESM2.xlsx",
                              sheet_name='NAS ordinal regression', engine='openpyxl')
gse130970_sig = gse130970_df[gse130970_df['adj_P'] < 0.01]
gse130970_symbols = set(gse130970_sig['gene_symbol'].dropna().astype(str))
gse130970_lfc = dict(zip(gse130970_sig['gene_symbol'].astype(str),
                          gse130970_sig['range_log2FC'].astype(float)))
gse130970_ensembl = symbols_to_ensembl(gse130970_symbols, ref_dvc['sym_map'])
datasets['GSE130970'] = {
    'label': 'GSE130970', 'n_samples': 78,
    'source': 'Published (ordinal)', 'contrast': 'Disease vs Control',
    'genes_ensembl': gse130970_ensembl,
    'n_study_degs': len(gse130970_ensembl),
    'lfc_by_symbol': gse130970_lfc,
    'ref': ref_dvc,
}
print(f"    {len(gse130970_symbols)} symbols, {len(gse130970_ensembl)} mapped")

print("  GSE135251 (Govaere) — per-study...")
gse135251_ps = pd.read_csv(f"{PER_STUDY}/GSE135251_de_results.csv")
gse135251_ps['gene_base'] = gse135251_ps['gene'].str.replace(r'\.\d+$', '', regex=True)
gse135251_sig = gse135251_ps[gse135251_ps['adj.P.Val'] < PERSTUDY_PADJ]
datasets['GSE135251'] = {
    'label': 'GSE135251', 'n_samples': 206,
    'source': 'Per-study (limma)', 'contrast': 'Disease vs Control',
    'genes_ensembl': set(gse135251_sig['gene_base']),
    'n_study_degs': len(gse135251_sig),
    'lfc_by_ensembl': dict(zip(gse135251_sig['gene_base'], gse135251_sig['logFC'])),
    'ref': ref_dvc,
}
print(f"    {len(gse135251_sig)} DEGs")

print("  GSE126848 (Suppli) — reproduced...")
gse126848 = pd.read_csv(f"{PUB_DEGS}/GSE126848/deseq2_NAFL_NASH_vs_Normalweight.csv", index_col=0)
gse126848.index = gse126848.index.astype(str).str.replace('"', '')
gse126848['gene_base'] = gse126848.index.str.replace(r'\.\d+$', '', regex=True)
gse126848_sig = gse126848[gse126848['padj'] < 0.05].copy()
datasets['GSE126848'] = {
    'label': 'GSE126848', 'n_samples': 31,
    'source': 'Reproduced (DESeq2)', 'contrast': 'Disease vs Control',
    'genes_ensembl': set(gse126848_sig['gene_base']),
    'n_study_degs': len(gse126848_sig),
    'lfc_by_ensembl': dict(zip(gse126848_sig['gene_base'], gse126848_sig['log2FoldChange'])),
    'ref': ref_dvc,
}
print(f"    {len(gse126848_sig)} DEGs")

# PRJNA512027 (Gerhard 2018) intentionally omitted from this comparison:
# L0/S0 library-prep / diagnosis confound — see top-of-file note.

print("  GSE213621 (Chen) — per-study...")
gse213621_ps = pd.read_csv(f"{PER_STUDY}/GSE213621_de_results.csv")
gse213621_ps['gene_base'] = gse213621_ps['gene'].str.replace(r'\.\d+$', '', regex=True)
gse213621_sig = gse213621_ps[gse213621_ps['adj.P.Val'] < PERSTUDY_PADJ]
datasets['GSE213621'] = {
    'label': 'GSE213621', 'n_samples': 368,
    'source': 'Per-study (limma)', 'contrast': 'Disease vs Control',
    'genes_ensembl': set(gse213621_sig['gene_base']),
    'n_study_degs': len(gse213621_sig),
    'lfc_by_ensembl': dict(zip(gse213621_sig['gene_base'], gse213621_sig['logFC'])),
    'ref': ref_dvc,
}
print(f"    {len(gse213621_sig)} DEGs")

print("  GSE162694 (Pantano) — per-study...")
gse162694_ps = pd.read_csv(f"{PER_STUDY}/GSE162694_de_results.csv")
gse162694_ps['gene_base'] = gse162694_ps['gene'].str.replace(r'\.\d+$', '', regex=True)
gse162694_sig = gse162694_ps[gse162694_ps['adj.P.Val'] < PERSTUDY_PADJ]
datasets['GSE162694'] = {
    'label': 'GSE162694', 'n_samples': 143,
    'source': 'Per-study (limma)', 'contrast': 'Disease vs Control',
    'genes_ensembl': set(gse162694_sig['gene_base']),
    'n_study_degs': len(gse162694_sig),
    'lfc_by_ensembl': dict(zip(gse162694_sig['gene_base'], gse162694_sig['logFC'])),
    'ref': ref_dvc,
}
print(f"    {len(gse162694_sig)} DEGs")

# --- NASH vs NAFL dataset (1) ---

print("  GSE167523 (Kodama) — NASH vs NAFL...")
gse167523 = pd.read_csv(f"{PUB_DEGS}/GSE167523/table_s2_genes.csv")
gse167523_symbols = set(gse167523['symbol'].dropna().astype(str))
gse167523_lfc = dict(zip(gse167523['symbol'].astype(str),
                         np.log2(gse167523['fold_change'].astype(float))))
gse167523_ensembl = symbols_to_ensembl(gse167523_symbols, ref_nn['sym_map'])
datasets['GSE167523'] = {
    'label': 'GSE167523', 'n_samples': 98,
    'source': 'Published', 'contrast': 'NASH vs NAFL',
    'genes_ensembl': gse167523_ensembl,
    'n_study_degs': len(gse167523_ensembl),
    'lfc_by_symbol': gse167523_lfc,
    'ref': ref_nn,
}
print(f"    {len(gse167523_symbols)} symbols, {len(gse167523_ensembl)} mapped")

# --- Adv vs Early Fibrosis datasets (3) ---
# Use script 15b per-study results (binary Advanced/Early, F2 dropped)
# These match the 15b mega/meta reference exactly

fib_ps = pd.read_csv(f"{DISEASE_SIG}/adv_vs_early_fibrosis_per_study.csv")

print("  GSE174478 (Kawamura) — adv vs early fib...")
gse174478_sub = fib_ps[fib_ps['dataset'] == 'GSE174478'].copy()
gse174478_sub['gene_base'] = gse174478_sub['gene'].str.replace(r'\.\d+$', '', regex=True)
gse174478_sig = gse174478_sub[gse174478_sub['adj.P.Val'] < PERSTUDY_PADJ]
datasets['GSE174478'] = {
    'label': 'GSE174478', 'n_samples': 94,
    'source': 'Per-study (limma)', 'contrast': 'Adv vs Early Fib',
    'genes_ensembl': set(gse174478_sig['gene_base']),
    'n_study_degs': len(gse174478_sig),
    'lfc_by_ensembl': dict(zip(gse174478_sig['gene_base'], gse174478_sig['logFC'])),
    'ref': ref_fib,
}
print(f"    {len(gse174478_sig)} DEGs")

print("  GSE193066 (Fujiwara) — adv vs early fib...")
gse193066_sub = fib_ps[fib_ps['dataset'] == 'GSE193066'].copy()
gse193066_sub['gene_base'] = gse193066_sub['gene'].str.replace(r'\.\d+$', '', regex=True)
gse193066_sig = gse193066_sub[gse193066_sub['adj.P.Val'] < PERSTUDY_PADJ]
datasets['GSE193066'] = {
    'label': 'GSE193066', 'n_samples': 164,
    'source': 'Per-study (limma)', 'contrast': 'Adv vs Early Fib',
    'genes_ensembl': set(gse193066_sig['gene_base']),
    'n_study_degs': len(gse193066_sig),
    'lfc_by_ensembl': dict(zip(gse193066_sig['gene_base'], gse193066_sig['logFC'])),
    'ref': ref_fib,
}
print(f"    {len(gse193066_sig)} DEGs")

print("  GSE240729 (Verschuren) — adv vs early fib...")
gse240729_sub = fib_ps[fib_ps['dataset'] == 'GSE240729'].copy()
gse240729_sub['gene_base'] = gse240729_sub['gene'].str.replace(r'\.\d+$', '', regex=True)
gse240729_sig = gse240729_sub[gse240729_sub['adj.P.Val'] < PERSTUDY_PADJ]
datasets['GSE240729'] = {
    'label': 'GSE240729', 'n_samples': 67,
    'source': 'Per-study (limma)', 'contrast': 'Adv vs Early Fib',
    'genes_ensembl': set(gse240729_sig['gene_base']),
    'n_study_degs': len(gse240729_sig),
    'lfc_by_ensembl': dict(zip(gse240729_sig['gene_base'], gse240729_sig['logFC'])),
    'ref': ref_fib,
}
print(f"    {len(gse240729_sig)} DEGs")

# ============================================================
# Compute overlaps using MATCHED references
# ============================================================
print("\nComputing overlaps (contrast-matched)...")

results = []
for ds_id, ds in datasets.items():
    genes = ds['genes_ensembl']
    ref = ds['ref']

    recovered_mega = genes & ref['mega_degs']
    recovered_meta = genes & ref['meta_degs']
    not_recovered = genes - (ref['mega_degs'] | ref['meta_degs'])
    meta_only = (genes & ref['meta_degs']) - ref['mega_degs']

    pct_mega = 100 * len(recovered_mega) / len(genes) if genes else 0
    pct_meta = 100 * len(recovered_meta) / len(genes) if genes else 0

    # LFC correlation
    lfc_pairs = []
    if 'lfc_by_symbol' in ds:
        for sym, pub_lfc in ds['lfc_by_symbol'].items():
            if sym in ref['sym_lfc'] and pd.notna(ref['sym_lfc'][sym]) and pd.notna(pub_lfc):
                lfc_pairs.append((float(pub_lfc), float(ref['sym_lfc'][sym])))
    elif 'lfc_by_ensembl' in ds:
        for ens, pub_lfc in ds['lfc_by_ensembl'].items():
            if ens in ref['mega_lfc'] and pd.notna(ref['mega_lfc'][ens]) and pd.notna(pub_lfc):
                lfc_pairs.append((float(pub_lfc), float(ref['mega_lfc'][ens])))

    if len(lfc_pairs) >= 10:
        pub_lfcs, mega_lfcs = zip(*lfc_pairs)
        r, _ = stats.pearsonr(pub_lfcs, mega_lfcs)
        concordant = sum(1 for pl, dl in lfc_pairs if np.sign(pl) == np.sign(dl))
        dir_conc = 100 * concordant / len(lfc_pairs)
    else:
        r, dir_conc = np.nan, np.nan
        pub_lfcs, mega_lfcs = [], []

    results.append({
        'dataset': ds_id, 'label': ds['label'], 'n_samples': ds['n_samples'],
        'source': ds['source'], 'contrast': ds['contrast'],
        'n_study_degs': ds['n_study_degs'],
        'n_recovered_mega': len(recovered_mega),
        'n_recovered_meta': len(recovered_meta),
        'n_meta_only': len(meta_only),
        'n_not_recovered': len(not_recovered),
        'pct_mega': pct_mega, 'pct_meta': pct_meta,
        'lfc_r': r, 'direction_concordance': dir_conc,
        'pub_lfcs': list(pub_lfcs) if lfc_pairs else [],
        'mega_lfcs': list(mega_lfcs) if lfc_pairs else [],
    })
    r_str = f"{r:.3f}" if not np.isnan(r) else "N/A"
    d_str = f"{dir_conc:.0f}%" if not np.isnan(dir_conc) else "N/A"
    print(f"  {ds_id} ({ds['label']}) [{ds['contrast']}]: "
          f"{len(genes):,} DEGs → mega={pct_mega:.1f}%, meta={pct_meta:.1f}%, r={r_str}, dir={d_str}")

df = pd.DataFrame(results)

# ============================================================
# Sort by contrast group, then by n_samples within group
# ============================================================
contrast_order = {'Disease vs Control': 0, 'NASH vs NAFL': 1, 'Adv vs Early Fib': 2}
df['contrast_rank'] = df['contrast'].map(contrast_order)
df = df.sort_values(['contrast_rank', 'n_samples'], ascending=[True, True]).reset_index(drop=True)

# ============================================================
# PLOTTING — Publication theme (Nature-compatible)
# ============================================================
print("\nGenerating figure...")

mpl.rcParams.update({
    'pdf.fonttype': 42, 'ps.fonttype': 42,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 7,
    'axes.titlesize': 8, 'axes.labelsize': 7,
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'axes.linewidth': 0.3,
    'xtick.major.width': 0.3, 'ytick.major.width': 0.3,
    'xtick.major.size': 2, 'ytick.major.size': 2,
    'axes.spines.top': False, 'axes.spines.right': False,
})

# MASLD palette
COL_MEGA     = '#0D47A1'   # deep blue
COL_META     = '#7B1FA2'   # violet
COL_NOT      = '#BDBDBD'   # grey
COL_PUB      = '#1565C0'   # blue
COL_REPRO    = '#00695C'   # teal
COL_PERSTUDY = '#C2185B'   # magenta
COL_CONC     = '#0D47A1'
COL_DISC     = '#C2185B'

# Contrast group colors (for brackets/shading)
CONTRAST_COLORS = {
    'Disease vs Control': '#0D47A1',
    'NASH vs NAFL': '#7B1FA2',
    'Adv vs Early Fib': '#00695C',
}
CONTRAST_BG = {
    'Disease vs Control': '#E3F2FD',
    'NASH vs NAFL': '#F3E5F5',
    'Adv vs Early Fib': '#E0F2F1',
}
CONTRAST_SHORT = {
    'Disease vs Control': 'Disease vs Control',
    'NASH vs NAFL': 'NASH vs NAFL',
    'Adv vs Early Fib': 'Adv. vs Early Fibrosis',
}

SOURCE_COLORS = {
    'Published': COL_PUB,
    'Published (ordinal)': COL_PUB,
    'Reproduced (DESeq2)': COL_REPRO,
    'Per-study (limma)': COL_PERSTUDY,
}

fig = plt.figure(figsize=(7.09, 6.2))
gs = gridspec.GridSpec(2, 2, hspace=0.45, wspace=0.35,
                       left=0.12, right=0.97, top=0.93, bottom=0.06,
                       height_ratios=[1.1, 1])

# ============================================================
# Panel A: Recovery Rate — grouped by contrast
# ============================================================
ax_a = fig.add_subplot(gs[0, 0])
y_pos = np.arange(len(df))
sizes = np.sqrt(df['n_study_degs'].values) * 3.5

# Add contrast group background shading
prev_contrast = None
group_start = 0
for i, (_, row) in enumerate(df.iterrows()):
    if row['contrast'] != prev_contrast:
        if prev_contrast is not None:
            # Shade previous group
            bg = CONTRAST_BG[prev_contrast]
            ax_a.axhspan(group_start - 0.45, i - 0.55, color=bg, alpha=0.4, zorder=0)
        group_start = i
        prev_contrast = row['contrast']
# Shade last group
ax_a.axhspan(group_start - 0.45, len(df) - 0.55, color=CONTRAST_BG[prev_contrast], alpha=0.4, zorder=0)

# Add horizontal separator lines between groups
group_boundaries = []
prev_contrast = None
for i, (_, row) in enumerate(df.iterrows()):
    if row['contrast'] != prev_contrast and prev_contrast is not None:
        group_boundaries.append(i - 0.5)
    prev_contrast = row['contrast']
for gb in group_boundaries:
    ax_a.axhline(gb, color='#9E9E9E', linewidth=0.4, linestyle='-', alpha=0.5, zorder=1)

for i, (_, row) in enumerate(df.iterrows()):
    c = SOURCE_COLORS.get(row['source'], '#999999')
    # Mega (large circle)
    ax_a.scatter(row['pct_mega'], i, s=sizes[i], c=c,
                 alpha=0.85, edgecolors='white', linewidths=0.4, zorder=3, marker='o')
    # Meta (diamond)
    ax_a.scatter(row['pct_meta'], i, s=20, c=c, marker='D',
                 alpha=0.7, edgecolors='white', linewidths=0.3, zorder=5)
    # Connecting line
    ax_a.plot([row['pct_meta'], row['pct_mega']], [i, i],
              color='#9E9E9E', linewidth=0.4, alpha=0.6, zorder=2)
    # DEG count annotation
    ax_a.text(max(row['pct_mega'], row['pct_meta']) + 2, i,
              f"{row['n_study_degs']:,}",
              va='center', fontsize=5, color='#616161')

# Contrast group labels on far right
group_ranges = {}
for i, (_, row) in enumerate(df.iterrows()):
    c = row['contrast']
    if c not in group_ranges:
        group_ranges[c] = [i, i]
    group_ranges[c][1] = i

# Median line
median_mega = df['pct_mega'].median()
ax_a.axvline(median_mega, color='#880E4F', linestyle='--', alpha=0.35, linewidth=0.5, zorder=1)
ax_a.text(median_mega + 1, len(df) - 0.3, f'median {median_mega:.0f}%',
          fontsize=5, color='#880E4F', alpha=0.7)

ax_a.set_yticks(y_pos)
ax_a.set_yticklabels([f"{row['label']}  (n={row['n_samples']})"
                       for _, row in df.iterrows()])
ax_a.set_xlabel('% Study DEGs recovered (padj < 0.1)')
ax_a.set_xlim(-3, 105)
ax_a.set_ylim(-0.6, len(df) - 0.4)
ax_a.set_title('a  Study DEG recovery by matched contrast', fontsize=8,
               fontweight='bold', loc='left', pad=4)

legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_PUB, markersize=4.5,
           markeredgecolor='white', markeredgewidth=0.3, label='Published'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_REPRO, markersize=4.5,
           markeredgecolor='white', markeredgewidth=0.3, label='Reproduced'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COL_PERSTUDY, markersize=4.5,
           markeredgecolor='white', markeredgewidth=0.3, label='Per-study (limma)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#757575', markersize=5,
           markeredgecolor='white', markeredgewidth=0.3, label='Mega (padj<0.1)'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='#757575', markersize=3.5,
           markeredgecolor='white', markeredgewidth=0.3, label='Meta (padj<0.1)'),
]
ax_a.legend(handles=legend_elements, loc='lower right', fontsize=4.5,
            frameon=True, framealpha=0.95, edgecolor='#E0E0E0',
            borderpad=0.4, handletextpad=0.3, ncol=1)

# ============================================================
# Panel B: LFC Concordance Scatter Grid — grouped by contrast
# ============================================================
gs_b = gridspec.GridSpecFromSubplotSpec(2, 5, subplot_spec=gs[0, 1],
                                         hspace=0.70, wspace=0.50)

for idx, (_, row) in enumerate(df.iterrows()):
    ax = fig.add_subplot(gs_b[idx // 5, idx % 5])
    pub_lfcs = np.array(row['pub_lfcs'])
    mega_lfcs = np.array(row['mega_lfcs'])
    contrast_col = CONTRAST_COLORS.get(row['contrast'], '#757575')

    if len(pub_lfcs) >= 10:
        concordant = np.sign(pub_lfcs) == np.sign(mega_lfcs)
        ax.scatter(pub_lfcs[concordant], mega_lfcs[concordant], s=0.3, alpha=0.15,
                   c=contrast_col, rasterized=True)
        ax.scatter(pub_lfcs[~concordant], mega_lfcs[~concordant], s=0.5, alpha=0.25,
                   c=COL_DISC, rasterized=True)
        lim = max(abs(pub_lfcs).max(), abs(mega_lfcs).max()) * 1.1
        lim = min(lim, 10)
        ax.plot([-lim, lim], [-lim, lim], color='#9E9E9E', linewidth=0.3, alpha=0.5)
        ax.axhline(0, color='#E0E0E0', linewidth=0.2)
        ax.axvline(0, color='#E0E0E0', linewidth=0.2)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)

        r_color = contrast_col if row['lfc_r'] > 0.5 else (COL_DISC if row['lfc_r'] < 0.3 else '#616161')
        ax.text(0.05, 0.95, f"r={row['lfc_r']:.2f}", transform=ax.transAxes,
                fontsize=5.5, fontweight='bold', va='top', color=r_color)
        ax.text(0.05, 0.75, f"{row['direction_concordance']:.0f}%",
                transform=ax.transAxes, fontsize=4.5, va='top', color='#757575')
    else:
        ax.text(0.5, 0.5, 'N/A', transform=ax.transAxes,
                ha='center', va='center', fontsize=6, color='#BDBDBD')

    # Title colored by contrast (accession label)
    title_name = row['label']
    ax.set_title(title_name, fontsize=5.5, fontweight='bold', pad=2,
                 color=contrast_col)
    ax.tick_params(labelsize=4, length=1.5, width=0.2)
    if idx % 5 == 0:
        ax.set_ylabel('Mega LFC', fontsize=5)
    if idx >= 5:
        ax.set_xlabel('Study LFC', fontsize=5)

fig.text(0.73, 0.945, 'b  LFC concordance (study vs matched mega)', fontsize=8,
         fontweight='bold', ha='center')

# ============================================================
# Panel C: Stacked Recovery Waterfall — grouped by contrast
# ============================================================
ax_c = fig.add_subplot(gs[1, 0])
# Sort by contrast group, then by n_study_degs within group
df_c = df.sort_values(['contrast_rank', 'n_study_degs'], ascending=[True, True]).reset_index(drop=True)
y = np.arange(len(df_c))
bar_height = 0.55

# Add contrast group background shading
prev_contrast = None
group_start = 0
for i, (_, row) in enumerate(df_c.iterrows()):
    if row['contrast'] != prev_contrast:
        if prev_contrast is not None:
            bg = CONTRAST_BG[prev_contrast]
            ax_c.axhspan(group_start - 0.45, i - 0.55, color=bg, alpha=0.4, zorder=0)
        group_start = i
        prev_contrast = row['contrast']
ax_c.axhspan(group_start - 0.45, len(df_c) - 0.55, color=CONTRAST_BG[prev_contrast], alpha=0.4, zorder=0)

# Separator lines
prev_contrast = None
for i, (_, row) in enumerate(df_c.iterrows()):
    if row['contrast'] != prev_contrast and prev_contrast is not None:
        ax_c.axhline(i - 0.5, color='#9E9E9E', linewidth=0.4, linestyle='-', alpha=0.5, zorder=1)
    prev_contrast = row['contrast']

bars_mega = df_c['n_recovered_mega'].values
bars_meta_only = df_c['n_meta_only'].values
bars_not = df_c['n_not_recovered'].values

ax_c.barh(y, bars_mega, height=bar_height, color=COL_MEGA, label='In mega', zorder=3)
ax_c.barh(y, bars_meta_only, height=bar_height, left=bars_mega,
          color=COL_META, alpha=0.7, label='In meta only', zorder=3)
ax_c.barh(y, bars_not, height=bar_height, left=bars_mega + bars_meta_only,
          color=COL_NOT, label='Not recovered', zorder=3)

for i, (_, row) in enumerate(df_c.iterrows()):
    total = row['n_study_degs']
    ax_c.text(total + total * 0.02, i, f"{row['pct_mega']:.0f}%",
              va='center', fontsize=5.5, fontweight='bold', color=COL_MEGA)

ax_c.set_yticks(y)
ax_c.set_yticklabels([row['label'] for _, row in df_c.iterrows()])
ax_c.set_xlabel('Number of DEGs')
ax_c.set_ylim(-0.6, len(df_c) - 0.4)
ax_c.set_title('c  Recovery counts (stacked)', fontsize=8, fontweight='bold', loc='left', pad=4)
ax_c.legend(loc='lower right', fontsize=5, frameon=True, framealpha=0.9,
            edgecolor='#E0E0E0', borderpad=0.4, handletextpad=0.3)

# ============================================================
# Panel D: Summary Heatmap — grouped by contrast
# ============================================================
ax_d = fig.add_subplot(gs[1, 1])

heatmap_cols = ['Study\nDEGs', '% Mega\nrecovered', '% Meta\nrecovered',
                'LFC r', 'Direction\nconcord.']
heatmap_data = []
contrast_labels = []
for _, row in df.iterrows():
    heatmap_data.append([
        row['n_study_degs'], row['pct_mega'], row['pct_meta'],
        row['lfc_r'] if not np.isnan(row['lfc_r']) else np.nan,
        row['direction_concordance'] if not np.isnan(row['direction_concordance']) else np.nan,
    ])
    contrast_labels.append(row['contrast'])

hm = np.array(heatmap_data, dtype=float)

# Normalize columns independently
hm_norm = np.zeros_like(hm)
for j in range(hm.shape[1]):
    col = hm[:, j]
    valid = ~np.isnan(col)
    if valid.sum() > 0:
        cmin, cmax = np.nanmin(col), np.nanmax(col)
        if cmax > cmin:
            hm_norm[:, j] = np.where(valid, (col - cmin) / (cmax - cmin), 0.5)
        else:
            hm_norm[:, j] = 0.5
hm_norm[np.isnan(hm)] = 0.5

cmap_blue = LinearSegmentedColormap.from_list('masld_blue',
    ['#FFFFFF', '#E3F2FD', '#90CAF9', '#42A5F5', '#1565C0', '#0D47A1'])

im = ax_d.imshow(hm_norm, cmap=cmap_blue, aspect='auto', vmin=0, vmax=1)

for i in range(hm.shape[0]):
    for j in range(hm.shape[1]):
        val = hm[i, j]
        if np.isnan(val):
            text = "N/A"
        elif j == 0:
            text = f"{int(val):,}"
        elif j in [1, 2, 4]:
            text = f"{val:.0f}%"
        else:
            text = f"{val:.2f}"
        bg = hm_norm[i, j]
        text_color = 'white' if bg > 0.65 else '#212121'
        ax_d.text(j, i, text, ha='center', va='center', fontsize=5.5,
                  fontweight='bold', color=text_color)

# Add contrast group brackets on the right
prev_contrast = None
group_ranges_d = {}
for i, cl in enumerate(contrast_labels):
    if cl not in group_ranges_d:
        group_ranges_d[cl] = [i, i]
    group_ranges_d[cl][1] = i

# Separator lines
prev_contrast = None
for i, cl in enumerate(contrast_labels):
    if cl != prev_contrast and prev_contrast is not None:
        ax_d.axhline(i - 0.5, color='#9E9E9E', linewidth=0.6, zorder=5)
    prev_contrast = cl

# Contrast bracket labels on right
for contrast_name, (start, end) in group_ranges_d.items():
    mid = (start + end) / 2
    short = {'Disease vs Control': 'DvC', 'NASH vs NAFL': 'NvN', 'Adv vs Early Fib': 'AvE'}[contrast_name]
    col = CONTRAST_COLORS[contrast_name]
    ax_d.text(hm.shape[1] - 0.5 + 0.8, mid, short, ha='center', va='center',
              fontsize=5.5, fontweight='bold', color=col,
              bbox=dict(boxstyle='round,pad=0.15', facecolor=CONTRAST_BG[contrast_name],
                        edgecolor=col, linewidth=0.3, alpha=0.8))

ax_d.set_xticks(np.arange(len(heatmap_cols)))
ax_d.set_xticklabels(heatmap_cols, fontsize=5.5)
ax_d.set_yticks(np.arange(len(df)))
ax_d.set_yticklabels([row['label'] for _, row in df.iterrows()])
ax_d.set_title('d  Summary (contrast-matched)', fontsize=8, fontweight='bold', loc='left', pad=4)
ax_d.tick_params(length=0)
ax_d.spines['left'].set_visible(False)
ax_d.spines['bottom'].set_visible(False)

# ============================================================
# Save
# ============================================================
os.makedirs(OUTDIR, exist_ok=True)

outpath_pdf = f"{OUTDIR}/published_deg_comparison.pdf"
fig.savefig(outpath_pdf, dpi=600, bbox_inches='tight', facecolor='white')
print(f"\nPDF saved: {outpath_pdf}")

summary_path = f"{OUTDIR}/published_deg_comparison_summary.csv"
df_out = df[['dataset', 'label', 'n_samples', 'source', 'contrast',
             'n_study_degs', 'n_recovered_mega', 'n_recovered_meta',
             'n_meta_only', 'n_not_recovered',
             'pct_mega', 'pct_meta', 'lfc_r', 'direction_concordance']].copy()
df_out.to_csv(summary_path, index=False)
print(f"Summary table: {summary_path}")

# Save LFC pairs for R figure integration
lfc_records = []
for ds_id, ds in datasets.items():
    ref = ds['ref']
    if 'lfc_by_symbol' in ds:
        for sym, pub_lfc in ds['lfc_by_symbol'].items():
            if sym in ref['sym_lfc'] and pd.notna(ref['sym_lfc'][sym]) and pd.notna(pub_lfc):
                lfc_records.append({
                    'dataset': ds_id, 'label': ds['label'],
                    'contrast': ds['contrast'],
                    'pub_lfc': float(pub_lfc), 'mega_lfc': float(ref['sym_lfc'][sym])
                })
    elif 'lfc_by_ensembl' in ds:
        for ens, pub_lfc in ds['lfc_by_ensembl'].items():
            if ens in ref['mega_lfc'] and pd.notna(ref['mega_lfc'][ens]) and pd.notna(pub_lfc):
                lfc_records.append({
                    'dataset': ds_id, 'label': ds['label'],
                    'contrast': ds['contrast'],
                    'pub_lfc': float(pub_lfc), 'mega_lfc': float(ref['mega_lfc'][ens])
                })
lfc_path = f"{OUTDIR}/published_deg_lfc_pairs.csv"
pd.DataFrame(lfc_records).to_csv(lfc_path, index=False)
print(f"LFC pairs: {lfc_path} ({len(lfc_records):,} pairs)")

plt.close()
print("\nDone!")
