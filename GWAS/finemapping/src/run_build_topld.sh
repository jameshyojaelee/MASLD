#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=220G
#SBATCH --time=48:00:00
#SBATCH --partition=bigmem
#SBATCH --array=0-87
#SBATCH --job-name=build_topld
#SBATCH --output=GWAS/finemapping/logs/build_topld_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/build_topld_%A_%a.err

# Build TOP-LD per-block LD matrices for all 4 ancestries × 22 chromosomes (88 tasks).
# Array index decomposition: ancestry = idx / 22; chr = (idx % 22) + 1.

set -o pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "${BASE}"

POPS=(EUR AFR EAS SAS)
IDX=${SLURM_ARRAY_TASK_ID}
POP_IDX=$((IDX / 22))
POP=${POPS[$POP_IDX]}
CHR=$((IDX % 22 + 1))

echo "[topld_build] task ${IDX}: ${POP} chr${CHR}  start $(date)"
Rscript GWAS/finemapping/src/build_topld_blocks.R "${POP}" "${CHR}"
echo "[topld_build] task ${IDX}: ${POP} chr${CHR}  done  $(date)"
