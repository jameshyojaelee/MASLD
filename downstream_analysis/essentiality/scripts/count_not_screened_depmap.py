import pandas as pd
import re

# File paths
depmap_path = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/essentiality/CRISPRGeneEffect.csv'
liang_path = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/essentiality/results/Essentiality_Quadrants_All.csv'

print("Loading DepMap data...")
# Read DepMap data
# It's large, but we only need column means.
# We can read in chunks or just read the whole thing if memory allows. 432MB is fine for pandas usually.
try:
    df_depmap = pd.read_csv(depmap_path, index_col=0) # Index is ModelID
except Exception as e:
    print(f"Error loading DepMap: {e}")
    exit(1)

print("Calculating mean scores...")
# Calculate mean essentiality score for each gene across all cell lines
gene_means = df_depmap.mean(axis=0)

# Parse gene names from "Gene (ID)" format
# Create a dictionary or series mapping clean names to scores
depmap_scores = {}
for col_name, score in gene_means.items():
    # Extract gene symbol (before the space and parenthesis)
    match = re.match(r"^(\S+)\s+\(\d+\)$", col_name)
    if match:
        gene_symbol = match.group(1)
        depmap_scores[gene_symbol] = score
    else:
        # Fallback if format is different or just gene name
        depmap_scores[col_name] = score

print(f"Total DepMap genes processed: {len(depmap_scores)}")

# Load Liang genes
print("Loading Liang screened genes...")
df_liang = pd.read_csv(liang_path)
# Genes are in 'Gene' column
liang_genes = set(df_liang['Gene'].tolist())
print(f"Total Liang genes: {len(liang_genes)}")

# Identify "Not Screened" genes
# These are genes in DepMap but NOT in Liang
not_screened_genes = {gene: score for gene, score in depmap_scores.items() if gene not in liang_genes}
print(f"Total 'Not Screened' genes: {len(not_screened_genes)}")

# Count genes with score between -1 and -0.5
# Range: -1.0 <= score <= -0.5
count_in_range = 0
genes_in_range = []

for gene, score in not_screened_genes.items():
    if -1.0 <= score <= -0.5:
        count_in_range += 1
        if len(genes_in_range) < 10: # Store first 10 for preview
            genes_in_range.append(f"{gene}: {score:.4f}")

print("-" * 30)
print(f"FINAL COUNT of 'Not Screened' genes with score [-1.0, -0.5]: {count_in_range}")
print("-" * 30)
print("Preview of genes:", genes_in_range)
