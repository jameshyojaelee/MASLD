import pandas as pd
import os
import sys

# Config
GSE_ID = "GSE159911"
ROOT_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"
METADATA_DIR = os.path.join(ROOT_DIR, GSE_ID, "metadata")
INPUT_FILE = os.path.join(METADATA_DIR, f"{GSE_ID}_SRARunTable.csv")
RUN_INFO_FILE = os.path.join(METADATA_DIR, f"{GSE_ID}_run_info.csv")
ATTR_FILE = os.path.join(METADATA_DIR, f"{GSE_ID}_sample_attributes.csv")

def main():
    if not os.path.exists(INPUT_FILE):
        print(f"Error: {INPUT_FILE} not found.")
        sys.exit(1)
        
    print(f"Reading {INPUT_FILE}...")
    df = pd.read_csv(INPUT_FILE)
    
    # 1. Create Run Info (Run, BioSample)
    # Pipeline expects columns: 'Run', 'BioSample' (case insensitive) -> handled by generate_samples_sheet.py
    # We'll keep standard columns
    if 'Run' not in df.columns or 'BioSample' not in df.columns:
        print("Error: Missing Run or BioSample columns.")
        sys.exit(1)
        
    run_df = df[['Run', 'BioSample']].copy()
    run_df.to_csv(RUN_INFO_FILE, index=False)
    print(f"Created {RUN_INFO_FILE}")
    
    # 2. Create Attributes (BioSample, Diet, Sex)
    # Pipeline expects 'BioSample' and columns like 'diet', 'treatment', 'sex'
    # We have 'Diet' in input.
    cols = ['BioSample', 'Diet']
    attr_df = df[cols].copy()
    
    # Add standardized columns if missing
    if 'Sex' not in df.columns:
        attr_df['Sex'] = 'Unspecified'
    else:
        attr_df['Sex'] = df['Sex']
        
    # Lowercase columns for pipeline consistency
    attr_df.columns = [c.lower() for c in attr_df.columns]
        
    attr_df.to_csv(ATTR_FILE, index=False)
    print(f"Created {ATTR_FILE}")

if __name__ == "__main__":
    main()
