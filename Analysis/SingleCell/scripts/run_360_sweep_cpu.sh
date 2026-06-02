#!/bin/bash
#SBATCH --job-name=s3_360_sweep
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --qos=nslab
#SBATCH --output=Analysis/SingleCell/scripts/logs/resolution_sweep/360_sweep_cpu_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/logs/resolution_sweep/360_sweep_cpu_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
mkdir -p Analysis/SingleCell/scripts/logs/resolution_sweep
export FORCE_CPU_LEIDEN=1
micromamba run -n rapids_singlecell python Analysis/SingleCell/scripts/360_resolution_sweep.py
