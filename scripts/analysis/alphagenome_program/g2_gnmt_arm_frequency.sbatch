#!/bin/bash
#SBATCH --job-name=agp_g-corrections
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_g-corrections-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_g-corrections-%j.err

# Recompute the four GNMT two-locus haplotype counts from data/1kg_eur by three routes, to settle
# the F1 (2/758) versus F2 (0.0000) disagreement. Reads only; writes one table into the
# g-integration package.

set -euo pipefail

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

export PYTHONNOUSERSITE=1
export G2_OUT_ROOT=$PROJ/GWAS/finemapping/results/alphagenome_program/g-integration-20260915T072500Z
export G2_WORK=$G2_OUT_ROOT/raw/gnmt_arm_frequency
export OMP_NUM_THREADS=1

cd "$PROJ"
echo "host $(hostname)  job ${SLURM_JOB_ID:-none}  start $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "cpus ${SLURM_CPUS_PER_TASK:-?}  mem ${SLURM_MEM_PER_NODE:-?}MB"

"$PY" "$PROJ/scripts/analysis/alphagenome_program/g2_gnmt_arm_frequency_resolution.py"

echo "end $(date -u +%Y-%m-%dT%H:%M:%SZ)"
