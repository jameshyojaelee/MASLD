
import pandas as pd
import re
import os

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
METADATA_DIR = os.path.join(BASE_DIR, "metadata")
SERIES_MATRIX = os.path.join(METADATA_DIR, "GSE167523_series_matrix.txt")
ENA_METADATA = os.path.join(METADATA_DIR, "full_metadata.txt")
OUTPUT_FILE = os.path.join(METADATA_DIR, "samples.tsv")

def parse_series_matrix(filepath):
    print(f"Parsing Series Matrix: {filepath}")
    metadata = {} 
    gsm_ids = []
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if line.startswith("!Sample_geo_accession"):
                parts = line.split('\t')
                gsm_ids = [x.strip('"') for x in parts[1:]]
                for gsm in gsm_ids:
                    metadata[gsm] = {}
            
            elif line.startswith("!Sample_title"):
                parts = line.split('\t')
                titles = [x.strip('"') for x in parts[1:]]
                for i, title in enumerate(titles):
                    if i < len(gsm_ids):
                        metadata[gsm_ids[i]]['title'] = title
                        
            elif line.startswith("!Sample_characteristics_ch1"):
                parts = line.split('\t')
                chars = [x.strip('"') for x in parts[1:]]
                for i, char in enumerate(chars):
                    if i < len(gsm_ids) and char:
                        if ':' in char:
                            key, val = char.split(':', 1)
                            # Standardize keys
                            key = key.strip().lower()
                            val = val.strip()
                            metadata[gsm_ids[i]][key] = val
                            
    return metadata

def parse_ena_metadata(filepath):
    print(f"Parsing ENA Metadata: {filepath}")
    mapping = {} 
    try:
        df = pd.read_csv(filepath, sep='\t')
        for idx, row in df.iterrows():
            srr = row['run_accession']
            gsm = None
            if 'sample_alias' in row and pd.notna(row['sample_alias']):
                match = re.search(r'(GSM\d+)', str(row['sample_alias']))
                if match: gsm = match.group(1)
            if not gsm and 'experiment_title' in row and pd.notna(row['experiment_title']):
                match = re.search(r'(GSM\d+)', str(row['experiment_title']))
                if match: gsm = match.group(1)
            if gsm:
                mapping[srr] = gsm
    except Exception as e:
        print(f"Error parsing ENA metadata: {e}")
    return mapping

def main():
    gsm_clinical = parse_series_matrix(SERIES_MATRIX)
    srr_to_gsm = parse_ena_metadata(ENA_METADATA)
    
    final_rows = []
    for srr, gsm in srr_to_gsm.items():
        if gsm in gsm_clinical:
            data = gsm_clinical[gsm]
            row = {
                'sample_id': srr,
                'gsm_id': gsm,
                'sample_name': data.get('title', ''),
                # Exact verification keys:
                # "disease state: NAFLD", "disease subtype: NAFL"
                # "age: 50", "gender: M"
                'condition': data.get('disease state', ''), 
                'subtype': data.get('disease subtype', ''),
                'sex': data.get('gender', ''),
                'age': data.get('age', ''),
                'tissue': data.get('tissue', ''),
            }
            final_rows.append(row)
            
    if final_rows:
        df_out = pd.DataFrame(final_rows)
        cols = ['sample_id', 'gsm_id', 'condition', 'subtype', 'age', 'sex', 'tissue', 'sample_name']
        df_out = df_out[cols]
        df_out.to_csv(OUTPUT_FILE, sep='\t', index=False)
        print(f"Saved standardized metadata to {OUTPUT_FILE}")
        print(df_out.head())
    else:
        print("Error: No merged data produced.")

if __name__ == "__main__":
    main()
