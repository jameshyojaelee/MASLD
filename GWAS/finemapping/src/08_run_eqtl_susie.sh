#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu
#SBATCH --array=1-22
#SBATCH --job-name=eqtl_susie
#SBATCH --output=logs/eqtl_susie_%A_%a.out
#SBATCH --error=logs/eqtl_susie_%A_%a.err

set -eo pipefail

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"

mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

# Pin the LD panel explicitly.  Without this the fits used whatever
# finemapping_functions.R defaulted to at run time, which is how the
# 2026-04 fits ended up on a retired layout with no record of it.
export LD_PANEL="${LD_PANEL:-polyfun}"
echo "LD_PANEL=${LD_PANEL}  EQTL_SUSIE_DIR=${EQTL_SUSIE_DIR:-<default>}"
export CHR=${SLURM_ARRAY_TASK_ID}
echo "Running eQTL SuSiE for chr${CHR}"
echo "Start: $(date)"

Rscript src/08_run_eqtl_susie.R

echo "Done: $(date)"
