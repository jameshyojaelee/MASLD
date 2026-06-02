#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=susiextopld
#SBATCH --output=GWAS/finemapping/logs/susiextopld_%A_%a.out
#SBATCH --error=GWAS/finemapping/logs/susiextopld_%A_%a.err

# Re-run SuSiEX joint EUR+EAS fine-mapping with LD_PANEL=topld so the EUR arm uses
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
export LD_PANEL=topld
export SUSIEX_RESULTS_SUFFIX="_topld"
# EAS arm should NOT swap to topld_eas (which lacks chr-level PLINK files);
# pin EAS to the 1kg_eas panel that has the .bed/.bim/.fam SuSiEX needs.
export EAS_LD_DIR="${FM_DIR}/data/ld_ref/1kg_eas"

echo "[susiextopld] LOCUS_ROW=${LOCUS_ROW}  LD_PANEL=${LD_PANEL}  start $(date)"
python src/10_run_susiex.py --locus-row "${LOCUS_ROW}" --threads 4
echo "[susiextopld] LOCUS_ROW=${LOCUS_ROW}  done  $(date)"
