#!/bin/bash
#SBATCH --job-name=s3_360_sweep
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/SingleCell/scripts/logs/resolution_sweep/360_sweep_cpu_clean_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/logs/resolution_sweep/360_sweep_cpu_clean_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
mkdir -p Analysis/SingleCell/scripts/logs/resolution_sweep
export FORCE_CPU_LEIDEN=1
# Avoid ~/.local fallback (has old anndata that can't read newer h5ad null encoding)
export PYTHONNOUSERSITE=1
# Use spatial env: scanpy 1.12 + sklearn 1.8 + anndata 0.12.10 (handles null encoding) + leidenalg 0.11
micromamba run -n spatial python -E -s Analysis/SingleCell/scripts/360_resolution_sweep.py
