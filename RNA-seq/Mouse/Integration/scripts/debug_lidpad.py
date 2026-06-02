import os
import csv
import gzip
import pandas as pd

BASE_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
gse159911_dir = os.path.join(BASE_DIR, "Mouse/Public_Diet_Models/GSE159911/metadata")
run_table_path = os.path.join(gse159911_dir, "GSE159911_SRARunTable.csv")
series_matrix_path = os.path.join(gse159911_dir, "GSE159911_series_matrix.txt.gz")

print(f"Checking paths:\n{run_table_path}\n{series_matrix_path}")

# 1. Map SRR (Run) -> GSM (Sample Name)
srr_to_gsm = {}
if os.path.exists(run_table_path):
    with open(run_table_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'Run' in row and 'Sample Name' in row:
                srr_to_gsm[row['Run']] = row['Sample Name']
    print(f"Loaded {len(srr_to_gsm)} SRR->GSM mappings.")
    print("Example SRR mapping:", list(srr_to_gsm.items())[:5])
else:
    print("Run Table not found!")

# 2. Map GSM -> Timepoint (from Series Matrix)
gsm_to_time = {}
if os.path.exists(series_matrix_path):
    with gzip.open(series_matrix_path, 'rt') as f:
        lines = f.readlines()
        
    geo_accessions = []
    time_characteristics = []
    
    for line in lines:
        line = line.strip()
        if line.startswith('!Sample_geo_accession'):
            geo_accessions = line.split('\t')[1:]
            geo_accessions = [g.strip('"') for g in geo_accessions]
        
        if line.startswith('!Sample_characteristics_ch1') and 'diet duration (weeks)' in line:
            vals = line.split('\t')[1:]
            time_characteristics = vals
            print(f"Captured time characteristics: {time_characteristics[:3]}...")
                
    if geo_accessions and time_characteristics:
         for i, gsm in enumerate(geo_accessions):
             if i < len(time_characteristics):
                 desc = time_characteristics[i].strip('"')
                 # desc like "diet duration (weeks): 8"
                 if 'diet duration (weeks):' in desc:
                     try:
                         val = desc.split(':')[-1].strip()
                         # Clean any extra text
                         val = val.split()[0] 
                         gsm_to_time[gsm] = f'{val}wk'
                     except:
                         pass
    
    print(f"Loaded {len(gsm_to_time)} GSM->Time mappings.")
    print("Example GSM mapping:", list(gsm_to_time.items())[:5])
    
    # Check intersection
    intersect = 0
    for srr, gsm in srr_to_gsm.items():
        if gsm in gsm_to_time:
            intersect += 1
            
    print(f"Intersection count (SRR -> GSM -> Time): {intersect}")
else:
    print("Series Matrix not found!")
