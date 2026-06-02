import pandas as pd
import re
import os

# Files
CORE_DEGS_FILE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/final_core_degs.csv"
DEPMAP_FILE = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/essentiality/CRISPRGeneEffect.csv'
LIANG_FILE = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/essentiality/results/Essentiality_Quadrants_All.csv'

print(f"Loading Core DEGs from {CORE_DEGS_FILE}...")
core_df = pd.read_csv(CORE_DEGS_FILE)

# Extract Human Symbols (Foreground)
# Need to handle 'human_ortholog_symbols' column standard usage seen in analyze_essentiality.py
foreground_genes = set()

# Simplified logic from analyze_essentiality.py
df_human = core_df.dropna(subset=['human_ortholog_symbols'])
df_human = df_human[df_human['human_ortholog_symbols'] != ""]

for _, row in df_human.iterrows():
    h_syms = str(row['human_ortholog_symbols']).split(';')
    for hs in h_syms:
        if hs.strip():
            foreground_genes.add(hs.strip())

print(f"Total Foreground Genes (Core Union): {len(foreground_genes)}")

# Load Liang Screened Genes
print(f"Loading Liang Data form {LIANG_FILE}...")
liang_df = pd.read_csv(LIANG_FILE)
liang_screened = set(liang_df['Gene'].tolist())
print(f"Liang Screened Genes: {len(liang_screened)}")

# Define 'Not Screened' in Foreground
not_screened_foreground = foreground_genes - liang_screened
print(f"Genes in Core Union but NOT in Liang: {len(not_screened_foreground)}")

# Load DepMap Data (Mean Scores)
print("Loading DepMap Data...")
# We need to map Gene Symbols to Scores
# Optimization: Read only columns that match our not_screened genes if possible, 
# but columns have IDs (Symbol (ID)). So typically need to read header first.
# Or just read all means like before (it worked fast enough).

df_depmap = pd.read_csv(DEPMAP_FILE, index_col=0)
gene_means = df_depmap.mean(axis=0)

depmap_scores = {}
for col_name, score in gene_means.items():
    # Extract gene symbol (before the space and parenthesis)
    match = re.match(r"^(\S+)\s+\(\d+\)$", col_name)
    if match:
        gene_symbol = match.group(1)
        depmap_scores[gene_symbol] = score
    else:
        depmap_scores[col_name] = score

# Get scores for 'Not Screened' genes
not_screened_scores = []
for gene in not_screened_foreground:
    if gene in depmap_scores:
        score = depmap_scores[gene]
        not_screened_scores.append({'Gene': gene, 'Score': score})
    else:
        # Gene not in DepMap either
        pass

df_not_screened = pd.DataFrame(not_screened_scores)
print(f"Found DepMap scores for {len(df_not_screened)} Not Screened genes.")

# Count Essential (< -0.5)
essential_all = df_not_screened[df_not_screened['Score'] < -0.5]
print(f"Probing 'Essential' (< -0.5): {len(essential_all)} (User saw 323?)")

# Count Target Range (-1.0 <= Score <= -0.5)
target_range = df_not_screened[(df_not_screened['Score'] >= -1.0) & (df_not_screened['Score'] <= -0.5)]
print(f"FINAL COUNT (Score [-1.0, -0.5]): {len(target_range)}")

print("Examples in range:", target_range.head(10)['Gene'].tolist())
