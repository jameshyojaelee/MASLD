#!/usr/bin/env bash
#SBATCH --job-name=MPRA
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --array=0-1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/MPRA_merge_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/MPRA_merge_%A_%a.err
set -euo pipefail
ROOT=${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}
cells=(HepG2 LX2)
export MASLD_PROJECT_ROOT=${ROOT}
python3 "${ROOT}/GWAS/finemapping/src/seqfunc_v2/mpra/01_merge_barcode_counts.py" \
  --cell-line "${cells[$SLURM_ARRAY_TASK_ID]}" --threads "${SLURM_CPUS_PER_TASK}"
