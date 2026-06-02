#!/bin/bash
#SBATCH --job-name=liver_netprox_lincs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=30:00:00
#SBATCH --output=logs/netprox_lincs_rerun_%j.out
#SBATCH --error=logs/netprox_lincs_rerun_%j.err

# Re-run Script 32 (network proximity, z-score fix) then Script 31 (LINCS annotation, filename fix)

set -euo pipefail

echo "=== Network Proximity + LINCS Re-run ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo ""

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/logs

echo "--- Running Script 32 (network proximity, z-score fix) ---"
Rscript RNA-seq/32_network_proximity.R

echo ""
echo "--- Running Script 31 (LINCS annotation, filename fix) ---"
Rscript RNA-seq/31_lincs_annotation.R

echo ""
echo "=== Complete ==="
echo "End: $(date)"
