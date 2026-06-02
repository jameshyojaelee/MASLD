import csv
import os
from pathlib import Path
import yaml

# Paths
ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq")
SAMPLESHEET = ROOT / "metadata/samples.tsv"
CONFIG = ROOT / "workflow/config.yaml"
# Note: config.yaml says "alignments/salmon", relative to project root
SALMON_OUT = ROOT / "alignments/salmon"

def verify_pipeline():
    # Load config to get excluded samples
    with open(CONFIG) as f:
        config_data = yaml.safe_load(f)
    
    excluded = set(config_data.get("excluded_samples", []))
    print(f"Loaded {len(excluded)} excluded samples from config.")

    # Load samples from metadata
    samples = []
    with open(SAMPLESHEET, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            samples.append(row['sample_id'])
    
    print(f"Total samples in metadata: {len(samples)}")
    
    # Check status
    missing_dir = []
    missing_quant = []
    completed = []
    skipped = []

    for sid in samples:
        if sid in excluded:
            skipped.append(sid)
            continue
            
        sample_dir = SALMON_OUT / sid
        quant_file = sample_dir / "quant.sf"
        
        if not sample_dir.exists():
            missing_dir.append(sid)
        elif not quant_file.exists():
            missing_quant.append(sid)
        else:
            completed.append(sid)
            
    print(f"Skipped (excluded): {len(skipped)}")
    print(f"Completed: {len(completed)}")
    print(f"Missing Directorires: {len(missing_dir)}")
    print(f"Missing quant.sf: {len(missing_quant)}")
    
    if missing_dir:
        print("\nSamples missing directory:")
        print(missing_dir[:20]) # Print first 20
    
    if missing_quant:
        print("\nSamples missing quant.sf:")
        print(missing_quant[:20])

if __name__ == "__main__":
    verify_pipeline()
