import os
import glob
import pandas as pd
import re

# Configuration
ROOT_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq"
OUTPUT_FILE = os.path.join(ROOT_DIR, "metadata", "samples.tsv")

def parse_sample_attributes(gse_dir):
    """
    Parses metadata files to map SRR -> {diet, sex}.
    Requires joining *_run_info.csv (SRR->BioSample) with *_sample_attributes.csv (BioSample->Attributes).
    """
    attr_file = glob.glob(os.path.join(gse_dir, "metadata", "*_sample_attributes.csv"))
    run_info_file = glob.glob(os.path.join(gse_dir, "metadata", "*_run_info.csv"))
    
    meta_map = {}
    
    if not attr_file or not run_info_file:
        return meta_map
    
    try:
        # Load Run Info (SRR -> BioSample)
        df_run = pd.read_csv(run_info_file[0])
        # Columns might be 'Run', 'BioSample'
        run_col = next((c for c in df_run.columns if c.lower() == 'run'), None)
        bio_col_run = next((c for c in df_run.columns if c.lower() == 'biosample'), None)
        
        if not run_col or not bio_col_run:
            print(f"Warning: Missing 'Run' or 'BioSample' in {run_info_file[0]}")
            return meta_map
            
        # Load Attributes (BioSample -> Diet, Sex)
        df_attr = pd.read_csv(attr_file[0])
        df_attr.columns = [c.lower() for c in df_attr.columns]
        
        # Identify potential condition columns
        bio_col_attr = next((c for c in df_attr.columns if 'biosample' in c or 'accession' in c), None)
        diet_col = next((c for c in df_attr.columns if 'diet' in c), None)
        treat_col = next((c for c in df_attr.columns if 'treatment' in c), None)
        sirna_col = next((c for c in df_attr.columns if 'sirna' in c), None)
        geno_col = next((c for c in df_attr.columns if 'genotype' in c), None)
        sex_col = next((c for c in df_attr.columns if 'sex' in c or 'gender' in c), None)
        
        if not bio_col_attr:
             print(f"Warning: Missing 'BioSample' in {attr_file[0]}")
             return meta_map
             
        # Merge
        # Map BioSample to Attributes
        attr_dict = {}
        for _, row in df_attr.iterrows():
            bs = str(row[bio_col_attr])
            
            # Construct Composite Condition
            parts = []
            if diet_col and str(row[diet_col]) not in ["nan", "None", "Unspecified"]:
                parts.append(str(row[diet_col]))
            if treat_col and str(row[treat_col]) not in ["nan", "None", "Unspecified"]:
                # Avoid redundancy if treatment == diet
                if not parts or parts[-1] != str(row[treat_col]):
                    parts.append(str(row[treat_col]))
            if sirna_col and str(row[sirna_col]) not in ["nan", "None", "Unspecified"]:
                parts.append(str(row[sirna_col]))
            # If empty, try genotype
            if not parts and geno_col:
                parts.append(str(row[geno_col]))
                
            condition = "_".join(parts) if parts else "Unspecified"
            sex = str(row[sex_col]) if sex_col else "Unspecified"
            
            # Clean
            condition = re.sub(r'[^a-zA-Z0-9_]', '_', condition)
            # Remove repeated underscores
            condition = re.sub(r'_+', '_', condition).strip('_')
            
            attr_dict[bs] = {'condition': condition, 'sex': sex}
            
        # Map SRR to Attributes via BioSample
        for _, row in df_run.iterrows():
            srr = str(row[run_col])
            bs = str(row[bio_col_run])
            
            if bs in attr_dict:
                meta_map[srr] = attr_dict[bs]
            else:
                 meta_map[srr] = {'condition': 'Unspecified', 'sex': 'Unspecified'}

    except Exception as e:
        print(f"Warning: Failed to parse metadata for {gse_dir}: {e}")
        
    return meta_map

def main():
    samples = []
    
    # Iterate over all GSE directories
    gse_dirs = glob.glob(os.path.join(ROOT_DIR, "GSE*"))
    
    for gse_dir in gse_dirs:
        gse_id = os.path.basename(gse_dir)
        fastq_dir = os.path.join(gse_dir, "fastq")
        
        if not os.path.isdir(fastq_dir):
            continue
            
        print(f"Processing {gse_id}...")
        
        # Get metadata mapping
        meta_map = parse_sample_attributes(gse_dir)
        
        # Find R1 files
        r1_files = glob.glob(os.path.join(fastq_dir, "*_1.fastq.gz"))
        
        for r1 in r1_files:
            # Infer R2
            r2 = r1.replace("_1.fastq.gz", "_2.fastq.gz")
            if not os.path.exists(r2):
                print(f"  Warning: Found {r1} but missing {r2}. Treating as single-end orphan?")
                continue
            
            # Extract SRR ID
            filename = os.path.basename(r1)
            srr_id = filename.split("_")[0]
            
            # Get metadata
            meta = meta_map.get(srr_id, {'condition': 'Unspecified', 'sex': 'Unspecified'})
            
            samples.append({
                'sample_id': srr_id,
                'fastq_1': r1,
                'fastq_2': r2,
                'strandedness': 'auto',
                'condition': meta['condition'],
                'sex': meta['sex'],
                'dataset': gse_id
            })

        # Logic for Single-End (files without _1/_2)
        # Exclude files that match _1.fastq.gz or _2.fastq.gz
        all_fastqs = glob.glob(os.path.join(fastq_dir, "*.fastq.gz"))
        single_end_files = [f for f in all_fastqs if "_1.fastq.gz" not in f and "_2.fastq.gz" not in f]

        for r1 in single_end_files:
            filename = os.path.basename(r1)
            # Assumption: filenames are SRRxxxx.fastq.gz
            if "_" in filename:
                 # fallback for safety, e.g. SRR_custom.fastq.gz
                 srr_id = filename.split("_")[0]
            else:
                 srr_id = filename.replace(".fastq.gz", "")

            meta = meta_map.get(srr_id, {'condition': 'Unspecified', 'sex': 'Unspecified'})
            
            samples.append({
                'sample_id': srr_id,
                'fastq_1': r1,
                'fastq_2': "", # Empty for single-end
                'strandedness': 'auto',
                'condition': meta['condition'],
                'sex': meta['sex'],
                'dataset': gse_id
            })
            
    # Convert to DataFrame and save
    if samples:
        df_samples = pd.DataFrame(samples)
        
        # Create metadata dir if explicit
        os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
        
        df_samples.to_csv(OUTPUT_FILE, sep='\t', index=False)
        print(f"Successfully generated {OUTPUT_FILE} with {len(df_samples)} samples.")
        
        # Preview
        print(df_samples.head())
        print(df_samples.groupby('dataset').size())
    else:
        print("No paired FASTQ files found.")

if __name__ == "__main__":
    main()
