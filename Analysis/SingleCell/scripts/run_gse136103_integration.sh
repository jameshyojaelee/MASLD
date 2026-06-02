#!/bin/bash
#SBATCH --job-name=gse136103_h5ad
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/gse136103_h5ad_%j.log

set -euo pipefail

echo "=== Phase 2: Convert GSE136103 processed matrices to h5ad ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

eval "$(micromamba shell hook --shell=bash)"
micromamba activate spatial

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p Analysis/SingleCell/scripts/logs

# Clean up any partial h5ad files from failed previous attempts
echo "Cleaning up partial GSM h5ad files..."
find Analysis/SingleCell/integration/input_h5ad/human/ -name 'GSM*.h5ad' -delete 2>/dev/null || true

python3 Analysis/SingleCell/integration/scripts/add_gse136103.py

echo ""
echo "=== Verifying output ==="
echo "H5ad files in input_h5ad/human/ matching GSM*:"
ls -lh Analysis/SingleCell/integration/input_h5ad/human/GSM*.h5ad 2>/dev/null | wc -l
echo ""
echo "Updated manifest entries for GSE136103:"
grep "GSE136103" Analysis/SingleCell/integration/input_h5ad/sample_manifest.csv | wc -l
echo ""
echo "Full manifest sample count:"
wc -l Analysis/SingleCell/integration/input_h5ad/sample_manifest.csv

echo ""
echo "End: $(date)"
