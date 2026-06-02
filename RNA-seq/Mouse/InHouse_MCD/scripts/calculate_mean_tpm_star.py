import pandas as pd
import os
import glob
from collections import defaultdict

# --- Configuration ---
BASE_DIR = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/in-house_MCD_RNAseq"
METADATA_FILE = os.path.join(BASE_DIR, "metadata/samples.tsv")
GTF_FILE = os.path.join(BASE_DIR, "reference/raw/gencode.vM33.annotation.gtf")
ALIGN_DIR = os.path.join(BASE_DIR, "alignments/star")
OUTPUT_DIR = BASE_DIR # Output to the main directory as implied by the request

def parse_gtf_gene_lengths(gtf_path):
    """
    Calculates the union exon length for each gene from a GTF file.
    Returns a dictionary: {gene_id: length} and {gene_id: gene_name}
    """
    print(f"Parsing GTF: {gtf_path}...")
    gene_exons = defaultdict(list)
    gene_names = {}
    
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            
            # Extract gene_id and gene_name using simple string parsing
            attributes = parts[8]
            gene_id = None
            gene_name = None
            
            # Helper to extract attribute value
            # Expecting format: gene_id "ENSMUSG..."; gene_name "Name";
            for attr in attributes.split(';'):
                attr = attr.strip()
                if not attr:
                    continue
                if ' ' in attr:
                    key, val = attr.split(' ', 1)
                    val = val.strip('"')
                    if key == 'gene_id':
                        gene_id = val
                    elif key == 'gene_name':
                        gene_name = val
            
            if not gene_id:
                continue
                
            if gene_name:
                gene_names[gene_id] = gene_name
            else:
                if gene_id not in gene_names:
                     gene_names[gene_id] = gene_id # Fallback
            
            if feature_type == 'exon':
                start = int(parts[3])
                end = int(parts[4])
                gene_exons[gene_id].append((start, end))

    print("Calculating union exon lengths...")
    gene_lengths = {}
    for gene_id, exons in gene_exons.items():
        # Merge overlapping intervals
        if not exons:
            gene_lengths[gene_id] = 0
            continue
            
        exons.sort()
        merged = []
        if exons:
            curr_start, curr_end = exons[0]
            for next_start, next_end in exons[1:]:
                if next_start < curr_end: # Overlap or adjacent (assuming 1-based closed, exact match < vs <= might differ but < is safe for overlap)
                    # Actually standard merge: if next_start <= curr_end + 1 (if we merge adjacent).
                    # But distinct exons usually imply distinct blocks. Let's merge strict overlaps: max(curr_end, next_end)
                    if next_start <= curr_end:
                         curr_end = max(curr_end, next_end)
                    else:
                        merged.append((curr_start, curr_end))
                        curr_start, curr_end = next_start, next_end
                else:
                    merged.append((curr_start, curr_end))
                    curr_start, curr_end = next_start, next_end
            merged.append((curr_start, curr_end))
        
        length = sum(end - start + 1 for start, end in merged)
        gene_lengths[gene_id] = length

    print(f"Processed {len(gene_lengths)} genes.")
    return gene_lengths, gene_names


def calculate_tpm(counts, lengths):
    """
    Calculates TPM from counts and lengths.
    counts: dict or series of {gene_id: count}
    lengths: dict of {gene_id: length_bp}
    """
    # Align indices
    common_genes = list(set(counts.keys()) & set(lengths.keys()))
    
    # Filter to common genes
    counts_vec = pd.Series({g: counts[g] for g in common_genes})
    lengths_vec = pd.Series({g: lengths[g] for g in common_genes})
    
    # RPK = Count / (Length / 1000)
    # Avoid division by zero by filtering length > 0 (which it should be)
    lengths_vec = lengths_vec[lengths_vec > 0]
    counts_vec = counts_vec[lengths_vec.index]
    
    rpk = counts_vec / (lengths_vec / 1000.0)
    
    # Scaling factor = Sum(RPK) / 1,000,000
    scaling_factor = rpk.sum() / 1000000.0
    
    if scaling_factor == 0:
        return pd.Series(0, index=counts_vec.index)
        
    tpm = rpk / scaling_factor
    return tpm

def main():
    # 1. Load Metadata
    print(f"Loading metadata from {METADATA_FILE}...")
    metadata = pd.read_csv(METADATA_FILE, sep='\t')
    # Expected columns: sample_id, week, diet, sex
    # Filter/Group
    control_samples = metadata[metadata['diet'] == 'Control']['sample_id'].tolist()
    mcd_samples = metadata[metadata['diet'] == 'MCD']['sample_id'].tolist()
    
    print(f"Control samples: {control_samples}")
    print(f"MCD samples: {mcd_samples}")
    
    # 2. Parse Gene Lengths
    gene_lengths, gene_names = parse_gtf_gene_lengths(GTF_FILE)
    
    # 3. Process All Samples
    all_tpms = {} # {sample_id: Series(tpm)}
    
    for sample_id in control_samples + mcd_samples:
        count_file = os.path.join(ALIGN_DIR, sample_id, f"{sample_id}.ReadsPerGene.out.tab")
        if not os.path.exists(count_file):
            print(f"WARNING: File not found for {sample_id}: {count_file}")
            continue
            
        print(f"Processing {sample_id}...")
        # STAR ReadsPerGene.out.tab format:
        # column 1: gene ID
        # column 2: counts for unstranded RNA-seq
        # column 3: counts for 1st read strand aligned with RNA (htseq-count option -s yes)
        # column 4: counts for 2nd read strand aligned with RNA (htseq-count option -s reverse)
        
        # We assume Column 4 (index 3) based on previous check (Reverse Stranded)
        df_counts = pd.read_csv(count_file, sep='\t', header=None, index_col=0, skiprows=4) 
        # skiprows=4 to skip N_unmapped, N_multimapping, N_noFeature, N_ambiguous
        
        counts_dict = df_counts.iloc[:, 2].to_dict() # Column index 2 is the 3rd data column (which is the 4th file column)
        
        tpm_series = calculate_tpm(counts_dict, gene_lengths)
        all_tpms[sample_id] = tpm_series

    # 4. Create DataFrame
    tpm_df = pd.DataFrame(all_tpms)
    # Add gene_names
    tpm_df['gene_name'] = tpm_df.index.map(lambda x: gene_names.get(x, x))
    
    # Move gene_name to front
    cols = ['gene_name'] + [c for c in tpm_df.columns if c != 'gene_name']
    tpm_df = tpm_df[cols]
    
    # 5. Calculate Means and Export
    
    # Control
    valid_controls = [s for s in control_samples if s in tpm_df.columns]
    if valid_controls:
        control_df = tpm_df[['gene_name'] + valid_controls].copy()
        control_df['mean_tpm'] = control_df[valid_controls].mean(axis=1)
        # Output columns: gene_id (index), gene_name, mean_tpm
        # We can just write the mean_tpm and gene_name if that's what is needed, 
        # or keeping it tidy with Index Name = gene_id
        control_out = control_df[['gene_name', 'mean_tpm']]
        control_out.index.name = 'gene_id'
        
        out_path = os.path.join(OUTPUT_DIR, "mean_tpm_control.csv")
        control_out.to_csv(out_path)
        print(f"Saved {out_path}")
    else:
        print("No valid control samples found.")

    # MCD
    valid_mcd = [s for s in mcd_samples if s in tpm_df.columns]
    if valid_mcd:
        mcd_df = tpm_df[['gene_name'] + valid_mcd].copy()
        mcd_df['mean_tpm'] = mcd_df[valid_mcd].mean(axis=1)
        
        mcd_out = mcd_df[['gene_name', 'mean_tpm']]
        mcd_out.index.name = 'gene_id'
        
        out_path = os.path.join(OUTPUT_DIR, "mean_tpm_mcd.csv")
        mcd_out.to_csv(out_path)
        print(f"Saved {out_path}")
    else:
        print("No valid MCD samples found.")
        
    print("Done.")

if __name__ == "__main__":
    main()
