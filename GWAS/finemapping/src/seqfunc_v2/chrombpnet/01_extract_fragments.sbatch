#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/extract_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/extract_%j.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT=${ROOT}
python3 "${ROOT}/GWAS/finemapping/src/seqfunc_v2/chrombpnet/01_extract_fragments.py"
python3 "${ROOT}/GWAS/finemapping/src/seqfunc_v2/chrombpnet/01b_audit_extraction.py"
