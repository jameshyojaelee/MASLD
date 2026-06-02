
import sys
import pandas as pd
from pathlib import Path
import os

def generate_samplesheet(metadata_file, fastq_dir, output_file):
    print(f"Reading metadata from: {metadata_file}")
    print(f"Looking for FASTQs in: {fastq_dir}")
    
    # Load verified metadata
    try:
        df = pd.read_csv(metadata_file, sep='\t')
    except Exception as e:
        print(f"Error reading metadata: {e}")
        sys.exit(1)
        
    fastq_path = Path(fastq_dir)
    if not fastq_path.exists():
        print(f"Error: FASTQ directory {fastq_dir} does not exist.")
        sys.exit(1)

    # Prepare new columns
    df['fastq_r1'] = ""
    df['fastq_r2'] = ""
    df['layout'] = ""

    # Iterate through samples and find files
    found_count = 0
    for idx, row in df.iterrows():
        sample_id = str(row['sample_id']).strip()
        
        # Search patterns
        # 1. Paired: {id}_1.fastq.gz / {id}_2.fastq.gz
        r1 = list(fastq_path.glob(f"{sample_id}_1.fastq.gz")) + list(fastq_path.glob(f"{sample_id}_1.fq.gz"))
        r2 = list(fastq_path.glob(f"{sample_id}_2.fastq.gz")) + list(fastq_path.glob(f"{sample_id}_2.fq.gz"))
        
        # 2. Single: {id}.fastq.gz
        single = list(fastq_path.glob(f"{sample_id}.fastq.gz")) + list(fastq_path.glob(f"{sample_id}.fq.gz"))
        
        if r1:
            df.at[idx, 'fastq_r1'] = str(r1[0].resolve())
            if r2:
                df.at[idx, 'fastq_r2'] = str(r2[0].resolve())
                df.at[idx, 'layout'] = 'PAIRED'
            else:
                df.at[idx, 'layout'] = 'SINGLE' # R1 only found, but marked paired pattern? unlikely but possible
            found_count += 1
        elif single:
            df.at[idx, 'fastq_r1'] = str(single[0].resolve())
            df.at[idx, 'layout'] = 'SINGLE'
            found_count += 1
        else:
            print(f"Warning: No FASTQ found for {sample_id}")

    # Filter to only rows with found FASTQs (fastq_r1 not empty)
    df_found = df[df['fastq_r1'] != ""]

    print(f"Matched FASTQs for {len(df_found)} / {len(df)} samples.")
    
    # Save merged
    df_found.to_csv(output_file, sep='\t', index=False)
    print(f"Saved merged samplesheet to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python generate_samplesheet.py <metadata_tsv> <fastq_dir> <output_tsv>")
        sys.exit(1)
    generate_samplesheet(sys.argv[1], sys.argv[2], sys.argv[3])
