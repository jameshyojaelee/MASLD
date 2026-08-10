#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/logs/extract_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2/logs/extract_%j.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
BASE=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/gse281367_hepatocyte_5fold_v2
export MASLD_PROJECT_ROOT=${ROOT}
python3 "${ROOT}/GWAS/finemapping/src/seqfunc_v2/chrombpnet/01_extract_fragments.py" \
  --out "${BASE}" --workers "${SLURM_CPUS_PER_TASK}" --expected-donors 12
python3 "${ROOT}/GWAS/finemapping/src/seqfunc_v2/chrombpnet/01b_audit_extraction.py" \
  --out "${BASE}" --expected-donors 12
