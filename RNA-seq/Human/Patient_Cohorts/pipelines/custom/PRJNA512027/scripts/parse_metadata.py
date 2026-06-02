
import pandas as pd
import os

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
METADATA_DIR = os.path.join(BASE_DIR, "metadata")
SRA_TABLE = os.path.join(METADATA_DIR, "SraRunTable.csv")
OUTPUT_FILE = os.path.join(METADATA_DIR, "samples.tsv")

def main():
    print(f"Parsing SraRunTable: {SRA_TABLE}")
    
    try:
        df = pd.read_csv(SRA_TABLE)
        
        # Columns: Run, AGE, sex, disease, disease_stage
        
        df_out = pd.DataFrame()
        df_out['sample_id'] = df['Run']
        df_out['age'] = df['AGE']
        df_out['sex'] = df['sex']
        df_out['condition'] = df['disease']
        
        if 'disease_stage' in df.columns:
            df_out['fibrosis_stage'] = df['disease_stage'] # Map stage to fibrosis column
            
        df_out['tissue'] = df['tissue']
       
        print(f"Processed {len(df_out)} samples.")
        print(df_out.head())
        
        df_out.to_csv(OUTPUT_FILE, sep='\t', index=False)
        print(f"Saved to {OUTPUT_FILE}")
        
    except Exception as e:
        print(f"Error parsing SraRunTable: {e}")

if __name__ == "__main__":
    main()
