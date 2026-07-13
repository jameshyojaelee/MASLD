import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import matplotlib as mpl

# Sanjana Lab Publication Typography Standards
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['font.sans-serif'] = ["Helvetica", "Arial", "DejaVu Sans"]
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['axes.titlesize'] = 6
mpl.rcParams['axes.labelsize'] = 6
mpl.rcParams['font.size'] = 6
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['xtick.labelsize'] = 6
mpl.rcParams['ytick.labelsize'] = 6

# Paths
INPUT_CSV = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/causal_inference/causal_inference_summary.csv"
OUTPUT_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/figures"
GENE_INFO = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/references/gencode.v49.annotation.gene_info.csv"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load Data
df = pd.read_csv(INPUT_CSV)

# Load Mapping
if os.path.exists(GENE_INFO):
    gene_map = pd.read_csv(GENE_INFO, names=['gene_id', 'gene_name', 'gene_type'], header=None)
    gene_map['ensembl_id_base'] = gene_map['gene_id'].str.split('.').str[0]
    map_dict = dict(zip(gene_map['ensembl_id_base'], gene_map['gene_name']))
    df['gene_symbol'] = df['gene'].map(map_dict)
    df['display_name'] = df['gene_symbol'].fillna(df['gene'])
else:
    df['display_name'] = df['gene']  # Fallback

df_out = df.copy()
# Bring display_name to front
cols = ['display_name'] + [c for c in df_out.columns if c != 'display_name']
df_out = df_out[cols]
df_out.to_csv(INPUT_CSV.replace('.csv', '_with_symbols.csv'), index=False)

df = df.dropna(subset=['twas_z', 'twas_fdr']) # Require at least TWAS results

# Colors (MR_COLOR retired 2026-04-22 — MR ditched from paper)
TWAS_COLOR = '#D81B60' # Human (Pink/Purple family)
SIG_COLOR = '#E53935'
NS_COLOR = '#BDBDBD'

# ==============================================================================
# Plot 1: TWAS Z-score Distribution (Density)
# ==============================================================================
plt.figure(figsize=(5, 4))
sns.histplot(data=df, x='twas_z', bins=50, color=TWAS_COLOR, alpha=0.6, kde=True)
plt.axvline(0, color='black', linestyle='--', linewidth=0.5)
print("[caption] TWAS Liver Z-score Distribution")
plt.xlabel('S-PrediXcan Z-score')
plt.ylabel('Gene Count')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'causal_1_twas_distribution.pdf'), dpi=600)
plt.close()

# ==============================================================================
# Plot 2: MR Forest Plot -- REMOVED 2026-04-22 (MR ditched from paper)
# Plot 3: TWAS vs MR Scatter -- REMOVED 2026-04-22
# Plot 4: TWAS/MR Hit Counts -- simplified to TWAS-only count
# ==============================================================================
plt.figure(figsize=(5, 4))
twas_sig = df['twas_fdr'] < 0.05
sns.barplot(x=['TWAS Significant\n(FDR < 0.05)'], y=[twas_sig.sum()], palette=[TWAS_COLOR])
print("[caption] Causal Inference Hit Counts (TWAS only; MR ditched 2026-04-22)")
plt.ylabel('Number of Genes')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'causal_4_hit_counts.pdf'), dpi=600)
plt.close()

# ==============================================================================
# Plot 5: "Manhattan"-style TWAS Volcano
# ==============================================================================
plt.figure(figsize=(6, 5))
plt.scatter(df['twas_z'], -np.log10(df['twas_p'] + 1e-300), color=NS_COLOR, s=5, alpha=0.5)

# Highlight top TWAS
sig_twas = df[df['twas_fdr'] < 0.05]
plt.scatter(sig_twas['twas_z'], -np.log10(sig_twas['twas_p']), color=TWAS_COLOR, s=15, alpha=0.8)

# Label top 3 TWAS per standard De-Clutter protocol
top_3 = sig_twas.sort_values('twas_p').head(3)
for _, row in top_3.iterrows():
    plt.text(row['twas_z'] + 0.1, -np.log10(row['twas_p']), row['display_name'], fontsize=6)

plt.axhline(-np.log10(max(sig_twas['twas_p'].max(), 1e-10)), color='black', linestyle='--', linewidth=0.5)
print("[caption] TWAS Liver Volcano Map")
plt.xlabel('TWAS Z-score')
plt.ylabel('-log10(p-value)')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'causal_5_twas_volcano.pdf'), dpi=600)
plt.close()

print("All python figures regenerated successfully with Gene Symbols.")
