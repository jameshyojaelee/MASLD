import pandas as pd
import gzip
import os

# Path to a sample data file
data_path = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/streamlit_deg_explorer/data/mcd_week_pooled_combined.tsv.gz"

print(f"Reading {data_path}...")
try:
    with gzip.open(data_path, 'rt') as f:
        df = pd.read_csv(f, sep='\t')
    
    print("Columns:", df.columns.tolist())
    if 'gene_type' in df.columns:
        print("\nBiotype counts in 'gene_type':")
        print(df['gene_type'].value_counts().head())
    elif 'biotype' in df.columns:
        print("\nBiotype counts in 'biotype':")
        print(df['biotype'].value_counts().head())
    else:
        print("\nNo direct 'biotype' or 'gene_type' column found.")
        # Check if we can infer from other columns?
        # Standard DESeq2 output might not have it unless added.
        
except Exception as e:
    print("Error:", e)
