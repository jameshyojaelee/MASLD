import pandas as pd
import re
import os

# Files
CORE_DEGS_FILE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/final_core_degs.csv"
DEPMAP_FILE = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/essentiality/CRISPRGeneEffect.csv'
LIANG_FILE = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/essentiality/results/Essentiality_Quadrants_All.csv'

print(f"Loading Core DEGs from {CORE_DEGS_FILE}...")
core_df = pd.read_csv(CORE_DEGS_FILE)

# Extract Human Symbols (Foreground aka Core Union)
foreground_genes = set()
df_human = core_df.dropna(subset=['human_ortholog_symbols'])
df_human = df_human[df_human['human_ortholog_symbols'] != ""]

for _, row in df_human.iterrows():
    h_syms = str(row['human_ortholog_symbols']).split(';')
    for hs in h_syms:
        if hs.strip():
            foreground_genes.add(hs.strip())

print(f"Total Core Union Genes: {len(foreground_genes)}")

# ---------------------------------------------------------
# Part 1: Screened in Liang (All Quadrants)
# ---------------------------------------------------------
print(f"\nLoading Liang Data from {LIANG_FILE}...")
liang_df = pd.read_csv(LIANG_FILE)

# Filter Liang genes that are in our Core Union (Foreground)
# (Ideally all Liang genes are in Core Union, but let's be strict)
liang_in_core = liang_df[liang_df['Gene'].isin(foreground_genes)].copy()
print(f"Liang Genes in Core Union: {len(liang_in_core)} (Full Liang N={len(liang_df)})")

# Count in Range [-1.0, -0.5] based on Liang Score
target_liang = liang_in_core[
    (liang_in_core['Mean_Essentiality_Score'] >= -1.0) & 
    (liang_in_core['Mean_Essentiality_Score'] <= -0.5)
]

print(f"Screened (Liang) in Range [-1.0, -0.5]: {len(target_liang)}")

# ---------------------------------------------------------
# Part 2: Not Screened (DepMap Only)
# ---------------------------------------------------------
liang_gene_set = set(liang_df['Gene'].tolist())
not_screened_genes = foreground_genes - liang_gene_set
print(f"\nNot Screened Genes (in Core Union): {len(not_screened_genes)}")

print("Loading DepMap Data...")
df_depmap = pd.read_csv(DEPMAP_FILE, index_col=0)
gene_means = df_depmap.mean(axis=0)

depmap_scores = {}
for col_name, score in gene_means.items():
    match = re.match(r"^(\S+)\s+\(\d+\)$", col_name)
    if match:
        gene_symbol = match.group(1)
        depmap_scores[gene_symbol] = score
    else:
        depmap_scores[col_name] = score

# Get scores for Not Screened
ns_scores = []
for gene in not_screened_genes:
    if gene in depmap_scores:
        ns_scores.append(depmap_scores[gene])

# Count in Range [-1.0, -0.5] based on DepMap Score
count_ns_target = sum(-1.0 <= s <= -0.5 for s in ns_scores)
print(f"Not Screened (DepMap) in Range [-1.0, -0.5]: {count_ns_target}")

# ---------------------------------------------------------
# Summary
# ---------------------------------------------------------
total_target = len(target_liang) + count_ns_target
print("\n" + "="*40)
print(f"TOTAL Core Genes with Score [-1.0, -0.5]: {total_target}")
print(f"  - From Screened (Liang Score): {len(target_liang)}")
print(f"  - From Not Screened (DepMap Score): {count_ns_target}")
print("="*40)
