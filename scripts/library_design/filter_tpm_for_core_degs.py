import pandas as pd
import os

# Paths
BASE_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CORE_DEGS_FILE = os.path.join(BASE_DIR, "final_core_degs.csv")
MCD_TPM_FILE = os.path.join(BASE_DIR, "RNA-seq/Mouse/InHouse_MCD/mean_tpm_mcd.csv")
OUTPUT_FILE = os.path.join(BASE_DIR, "RNA-seq/Mouse/InHouse_MCD/final_core_degs_mean_tpm_mcd.csv")

def main():
    print(f"Loading Core DEGs from {CORE_DEGS_FILE}...")
    core_degs = pd.read_csv(CORE_DEGS_FILE)
    
    print(f"Loading MCD TPMs from {MCD_TPM_FILE}...")
    mcd_tpm = pd.read_csv(MCD_TPM_FILE)
    
    # Check columns
    print("Core DEGs columns:", core_degs.columns.tolist())
    print("MCD TPM columns:", mcd_tpm.columns.tolist())
    
    # Process Match Key: Ensembl ID
    # Core DEGs: 'mouse_gene_id' (e.g., ENSMUSG00000052131)
    # MCD TPM: 'gene_id' (e.g., ENSMUSG00000052131.5) -> Need to strip version
    
    mcd_tpm['gene_id_clean'] = mcd_tpm['gene_id'].apply(lambda x: x.split('.')[0])
    
    # Filter/Merge
    # We want to keep all genes in core_degs, and find their TPM.
    # Note: Some core DEGs might not be in the TPM file (unlikely if they are from STAR/GTF, but possible).
    
    merged = core_degs.merge(mcd_tpm[['gene_id_clean', 'mean_tpm']], 
                             left_on='mouse_gene_id', 
                             right_on='gene_id_clean', 
                             how='left')
    
    # Check for missing values
    missing = merged[merged['mean_tpm'].isna()]
    if not missing.empty:
        print(f"WARNING: {len(missing)} genes from Core DEGs did not match validation TPMs.")
        print(missing[['mouse_gene_id', 'mouse_gene_symbol']].head())
    
    # Create final output
    # User might want just the values, but context usually implies strict tabular correlation.
    # I'll output relevant ID/Symbol columns and the calculated mean TPM.
    
    output_df = merged[['mouse_gene_id', 'mouse_gene_symbol', 'mean_tpm']].copy()
    output_df.rename(columns={'mean_tpm': 'in_house_mcd_mean_tpm'}, inplace=True)
    
    # Fill NA with 0 if appropriate, or leave as NA. Assuming 0 if not found is safer for filtering? 
    # Actually, if it's missing from STAR output it might mean zero counts or missing annotation.
    # Given we parsed the same GTF, it should match. If NA, let's keep it explicit for now or fill 0.
    # Let's fill 0 as likely non-detected.
    output_df['in_house_mcd_mean_tpm'] = output_df['in_house_mcd_mean_tpm'].fillna(0)
    
    print(f"Writing output to {OUTPUT_FILE}...")
    output_df.to_csv(OUTPUT_FILE, index=False)
    print("Done.")

if __name__ == "__main__":
    main()
