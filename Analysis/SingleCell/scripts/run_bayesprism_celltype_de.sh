#!/bin/bash
#SBATCH --job-name=bp_ct_de
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/bp_ct_de_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs/bp_ct_de_%j.err
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00

set -eo pipefail

mkdir -p /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/scripts/logs
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== bp_ct_de @ $(hostname) ==="
echo "Started: $(date)"
echo "Job ID:  ${SLURM_JOB_ID:-NA}"

Rscript Analysis/SingleCell/scripts/run_bayesprism_celltype_de.R

echo ""
echo "Running summary script..."
Rscript Analysis/SingleCell/scripts/summarize_bayesprism_celltype_de.R

echo "Done: $(date)"
