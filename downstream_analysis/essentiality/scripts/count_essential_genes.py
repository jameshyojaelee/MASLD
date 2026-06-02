import pandas as pd

# Load the data
file_path = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/downstream_analysis/essentiality/results/Essentiality_Quadrants_All.csv'
df = pd.read_csv(file_path)

# Filter for essential genes (excluding Non_Essential)
essential_genes = df[df['Essentiality_Quadrant'] != 'Non_Essential']

# Filter for score between -0.5 and -1
# Note: Score is typically negative, so we look for range [-1, -0.5]
target_genes = essential_genes[
    (essential_genes['Mean_Essentiality_Score'] >= -1.0) & 
    (essential_genes['Mean_Essentiality_Score'] <= -0.5)
]

print(f"Total essential genes: {len(essential_genes)}")
print(f"Essential genes with score between -0.5 and -1: {len(target_genes)}")
print("Genes:", target_genes['Gene'].tolist())
print("\nDetails:")
print(target_genes[['Gene', 'Mean_Essentiality_Score', 'Essentiality_Quadrant']].to_string(index=False))
