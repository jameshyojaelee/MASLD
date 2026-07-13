#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/gate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/gate_%j.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
BASE=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2
python3 "${ROOT}/GWAS/finemapping/src/seqfunc_v2/chrombpnet/06_update_gate.py" \
  --internal "${BASE}/fold0/internal_qc.json" \
  --external "${BASE}/currin_calibration/currin_external_calibration.json"
