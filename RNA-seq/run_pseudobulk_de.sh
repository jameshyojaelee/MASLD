#!/bin/bash
#SBATCH --job-name=liver_pseudobulk_de
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/pseudobulk_de_%j.out
#SBATCH --error=logs/pseudobulk_de_%j.err

# Script 36: Pseudobulk DE from scVI-integrated atlas
# BLOCKED: Run only after scVI integration completes (produces integrated_atlas.h5ad)

set -euo pipefail

echo "=== Pseudobulk DE ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/logs

# Check that scVI output exists
H5AD="Analysis/SingleCell/results/integrated_atlas.h5ad"
if [ ! -f "$H5AD" ]; then
  echo "ERROR: $H5AD not found. scVI integration must complete first."
  exit 1
fi

echo "Running Script 36..."
Rscript RNA-seq/36_pseudobulk_de.R

echo ""
echo "=== Complete ==="
echo "End: $(date)"
