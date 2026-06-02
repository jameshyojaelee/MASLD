#!/usr/bin/env bash
#SBATCH --job-name=BayesPrism-attribution
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=210G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/bayesprism_attribution_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/bayesprism_attribution_%j.err

# ---------------------------------------------------------------------------
# Stage 2: BayesPrism unified attribution (Script 22)
#
# Depends on Stage 1 (21_bayesprism_kallisto.sbatch) completing.
# Submit with: sbatch --dependency=afterok:<STAGE1_JOBID> 22_bayesprism_attribution_kallisto.sbatch
# ---------------------------------------------------------------------------

set -eo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

set -u

export R_LIBS_USER=""

echo "=============================="
echo "BayesPrism Attribution (Script 22)"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "=============================="

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Back up existing results
OUTDIR="RNA-seq/Human/Patient_Cohorts/analysis/integration/results/deconvolution/bayesprism"
for f in "${OUTDIR}"/*.csv; do
  if [[ -f "$f" ]]; then
    backup="${f%.csv}_star_backup.csv"
    if [[ ! -f "$backup" ]]; then
      cp "$f" "$backup"
      echo "Backed up: $(basename $f) -> $(basename $backup)"
    fi
  fi
done

Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/22_bayesprism_deconvolution.R

echo ""
echo "Finished: $(date)"
echo "Outputs:"
ls -lh "${OUTDIR}"/ 2>/dev/null
