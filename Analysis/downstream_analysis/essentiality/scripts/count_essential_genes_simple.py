import pandas as pd

# Load the data
file_path = '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/downstream_analysis/essentiality/results/Essentiality_Quadrants_All.csv'
df = pd.read_csv(file_path)

# Filter for essential genes (excluding Non_Essential)
essential_genes = df[df['Essentiality_Quadrant'] != 'Non_Essential']

# Filter for score between -0.5 and -1
target_genes = essential_genes[
    (essential_genes['Mean_Essentiality_Score'] >= -1.0) & 
    (essential_genes['Mean_Essentiality_Score'] <= -0.5)
]

print(f"FINAL COUNT: {len(target_genes)}")
