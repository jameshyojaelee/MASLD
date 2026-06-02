import pandas as pd
import os

# Define file paths
BASE_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/InHouse_MCD"
METADATA_FILE = os.path.join(BASE_DIR, "metadata/samples.tsv")
COUNTS_FILE = os.path.join(BASE_DIR, "normalized_counts_all_samples_salmon.csv")
# UPDATED: Using the combined pooled deseq2 results
DE_FILE = os.path.join(BASE_DIR, "results/combined_pooled/deseq2_results.tsv")
OUTPUT_FILE = os.path.join(BASE_DIR, "MCD_mean_tpm_filtered_combined.csv")

def main():
    print("Loading metadata...")
    metadata = pd.read_csv(METADATA_FILE, sep="\t")
    
    # Identify MCD samples
    # Looking for samples where diet is MCD
    mcd_samples = metadata[metadata['diet'] == 'MCD']['sample_id'].tolist()
    print(f"Indices of MCD samples identified: {mcd_samples}")

    print("Loading normalized counts...")
    counts_df = pd.read_csv(COUNTS_FILE, index_col=0) # Assuming Gene is index
    # Ensure all MCD samples are in counts columns
    missing_samples = [s for s in mcd_samples if s not in counts_df.columns]
    if missing_samples:
        raise ValueError(f"Missing samples in counts file: {missing_samples}")

    print(f"Loading DE statistics from {DE_FILE}...")
    de_df = pd.read_csv(DE_FILE, sep="\t", index_col=0) # Assuming Gene is index
    
    # Filter genes based on criteria
    # Criteria: padj < 0.1 AND log2FoldChange > 0.8
    # Using 'padj' and 'log2FoldChange' columns from the deseq2 results file
    
    # Ensure columns exist
    if 'padj' not in de_df.columns or 'log2FoldChange' not in de_df.columns:
        raise ValueError("Required DE columns 'padj' or 'log2FoldChange' not found.")

    filtered_genes_mask = (de_df['padj'] < 0.1) & (de_df['log2FoldChange'] > 0.8)
    filtered_genes = de_df[filtered_genes_mask].index.tolist()
    
    print(f"Total genes in DE file: {len(de_df)}")
    print(f"Genes passing filter (padj < 0.1, log2FC > 0.8): {len(filtered_genes)}")

    # Subset counts to filtered genes and MCD samples
    # Use intersection of index to handle any missing genes
    common_genes = list(set(filtered_genes).intersection(set(counts_df.index)))
    print(f"Genes found in counts file: {len(common_genes)}")
    
    mcd_counts_filtered = counts_df.loc[common_genes, mcd_samples]

    # Calculate Mean TPM
    mean_tpm = mcd_counts_filtered.mean(axis=1)
    
    # Create result DataFrame
    result_df = pd.DataFrame({'Mean_TPM': mean_tpm})
    result_df.index.name = 'Gene'
    
    # Sort by Mean TPM descending for readability
    result_df = result_df.sort_values(by='Mean_TPM', ascending=False)

    print(f"Saving results to {OUTPUT_FILE}...")
    result_df.to_csv(OUTPUT_FILE)
    print("Done.")

if __name__ == "__main__":
    main()
