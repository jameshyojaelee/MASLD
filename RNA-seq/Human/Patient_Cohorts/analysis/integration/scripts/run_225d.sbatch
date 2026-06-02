#!/bin/bash
#SBATCH --job-name=225d_ordinal
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/225d_%j.out
#SBATCH --error=logs/225d_%j.err
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
eval "$(micromamba shell hook --shell bash)" && micromamba activate spatial
python 225d_plasma_sweep_ordinal.py
