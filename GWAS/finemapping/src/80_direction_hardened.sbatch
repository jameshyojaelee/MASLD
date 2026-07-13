#!/bin/bash
#SBATCH --job-name=direction
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/80_direction_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/80_direction_%j.err
set -eo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$ROOT/GWAS/finemapping"
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
export PYTHONWARNINGS=ignore OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} \
       OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
PANELS="${1:-eqtl,caqtl,sqtl}"
NPERM="${2:-1000}"
NBOOT="${3:-2000}"
echo "[80] panels=$PANELS nperm=$NPERM nboot=$NBOOT"
python src/80_direction_hardened.py --panels "$PANELS" --nperm "$NPERM" --nboot "$NBOOT"
echo "[80] done"
