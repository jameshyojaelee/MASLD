#!/bin/bash
#SBATCH --job-name=225b_svm
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/225b_svm_%j.out
#SBATCH --error=logs/225b_svm_%j.err
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
eval "$(micromamba shell hook --shell bash)" && micromamba activate spatial
export MODEL_SUBSET="linear_svc,rbf_svm,nu_svc"
export OUT_SUFFIX="svm"
python 225b_plasma_sweep_boosting.py
