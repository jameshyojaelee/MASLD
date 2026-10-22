#!/bin/bash
#SBATCH --job-name=fig2_native_proposal
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=logs/fig2_native_proposal_%j.out
#SBATCH --error=logs/fig2_native_proposal_%j.err
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export MASLD_PROJECT_ROOT="$PWD"
export FIG2_PROPOSAL_ROOT="$PWD/figures/candidates/fig2-native-review-proposal-20261014_${SLURM_JOB_ID}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python scripts/figures/assemble_fig2_native_proposal.py
pdftotext "$FIG2_PROPOSAL_ROOT/Figure2_A-I_native_review_PROPOSAL.pdf" "$FIG2_PROPOSAL_ROOT/Figure2_A-I_native_review_PROPOSAL.txt"
pdftoppm -r 150 -singlefile -png "$FIG2_PROPOSAL_ROOT/Figure2_A-I_native_review_PROPOSAL.pdf" "$FIG2_PROPOSAL_ROOT/Figure2_A-I_native_review_PROPOSAL"
cp scripts/figures/assemble_fig2_native_proposal.py scripts/figures/run_fig2_native_proposal.sbatch "$FIG2_PROPOSAL_ROOT/"
