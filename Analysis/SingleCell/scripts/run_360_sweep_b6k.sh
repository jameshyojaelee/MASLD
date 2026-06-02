#!/bin/bash
#SBATCH --job-name=s3_360_sweep
#SBATCH --partition=gpu
#SBATCH --gres=gpu:b6k:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --qos=nslab
#SBATCH --output=Analysis/SingleCell/scripts/logs/resolution_sweep/360_sweep_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/logs/resolution_sweep/360_sweep_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
mkdir -p Analysis/SingleCell/scripts/logs/resolution_sweep
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH:-}:/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/lib/python3.12/site-packages/nvidia/cu13/lib
micromamba run -n rapids_singlecell python Analysis/SingleCell/scripts/360_resolution_sweep.py
