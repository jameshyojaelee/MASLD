#!/bin/bash
#SBATCH --job-name=stg162_loco_nmf
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/logs/162_loco_nmf_%j.out
#SBATCH --error=scripts/logs/162_loco_nmf_%j.err

source ~/.bashrc
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

python scripts/162_loco_nmf_subtypes.py
