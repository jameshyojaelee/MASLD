#!/bin/bash
#SBATCH --job-name=fc_merge_GSE240729
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/logs/fc_merge_%j.log

set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

OUTDIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE240729/counts/featurecounts"
PERDIR="${OUTDIR}/per_sample"

echo "Merging per-sample featureCounts..."

COUNT=$(ls "${PERDIR}"/*.counts.txt 2>/dev/null | wc -l)
echo "Found ${COUNT}/66 per-sample count files"
if [ "$COUNT" -lt 66 ]; then
    echo "ERROR: Missing count files. Expected 66, found ${COUNT}"
    exit 1
fi

python3 - "${PERDIR}" "${OUTDIR}/gene_counts.txt" << 'PYEOF'
import sys, os, glob
import pandas as pd
from pathlib import Path

per_dir = sys.argv[1]
out_file = sys.argv[2]

count_files = sorted(glob.glob(os.path.join(per_dir, "*.counts.txt")))
print(f"Merging {len(count_files)} files...")

first = pd.read_csv(count_files[0], sep="\t", comment="#")
annot_cols = first.columns[:6].tolist()
merged = first[annot_cols].copy()

for cf in count_files:
    df = pd.read_csv(cf, sep="\t", comment="#")
    sample_name = Path(cf).stem.replace(".counts", "")
    merged[sample_name] = df.iloc[:, -1]

with open(out_file, "w") as fh:
    fh.write(f"# Program:featureCounts v2.1.1; Merged from {len(count_files)} per-sample runs\n")
merged.to_csv(out_file, sep="\t", index=False, mode="a")

print(f"Output: {out_file}")
print(f"Dimensions: {merged.shape[0]} genes x {merged.shape[1] - 6} samples ({merged.shape[1]} total cols)")
PYEOF

echo "Done: $(date)"
