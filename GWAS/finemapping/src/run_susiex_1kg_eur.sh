#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=susiex1kg
#SBATCH --output=GWAS/finemapping/logs/susiex1kg_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/susiex1kg_%A_%a.err

# Re-run SuSiEX joint EUR+EAS fine-mapping with LD_PANEL=1kg so the EUR arm uses
# our new 1kg_eur panel instead of sghatan UKBB EUR. EAS arm continues to use
# 1kg_eas (unchanged — no sghatan dependency). Output routed to a _1kg-suffixed
# directory via RESULTS_SUFFIX (handled by the Python script below).

set -o pipefail
export PYTHONNOUSERSITE=1

FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs

eval "$(micromamba shell hook -s bash)"
micromamba activate susiex
set -u

LOCUS_ROW=${SLURM_ARRAY_TASK_ID:-1}
export LD_PANEL=1kg
export SUSIEX_RESULTS_SUFFIX="_1kg"

echo "[susiex1kg] LOCUS_ROW=${LOCUS_ROW}  LD_PANEL=${LD_PANEL}  start $(date)"
python src/10_run_susiex.py --locus-row "${LOCUS_ROW}" --threads 4
echo "[susiex1kg] LOCUS_ROW=${LOCUS_ROW}  done  $(date)"
