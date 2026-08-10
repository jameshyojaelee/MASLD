#!/bin/bash
#SBATCH --job-name=susieR
#SBATCH --partition=cpu
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/reorient_audit/run_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/tmp_diag/reorient_audit/run_%j.err
set +u
export PATH="${HOME}/.local/bin:${PATH}"
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
Rscript tmp_diag/reorient_audit/extract_and_consensus.R
