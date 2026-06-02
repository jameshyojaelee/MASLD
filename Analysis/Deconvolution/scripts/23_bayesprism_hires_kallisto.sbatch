#!/usr/bin/env bash
#SBATCH --job-name=BayesPrism-hires
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=210G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/bayesprism_hires_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/logs/bayesprism_hires_%j.err

# ---------------------------------------------------------------------------
# Stage 3: BayesPrism high-resolution cell-type expression (Script 118b)
#
# Runs on merged_dge.rds directly (847 samples, all 5 mega cohorts).
# Does NOT depend on per-cohort BayesPrism (Stage 1) -- runs its own
# full BayesPrism from scratch on the merged data.
#
# Depends on Stage 1 only logically (same kallisto counts era).
# Submit with: sbatch --dependency=afterok:<STAGE1_JOBID> 23_bayesprism_hires_kallisto.sbatch
# ---------------------------------------------------------------------------

set -eo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

set -u

export R_LIBS_USER=""

echo "=============================="
echo "BayesPrism High-Resolution (Script 118b)"
echo "Host: $(hostname)"
echo "Started: $(date)"
echo "=============================="

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

# Back up existing results
OUTDIR="RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/cibersortx_celltype_expression"
for f in "${OUTDIR}"/bayesprism_*.csv.gz "${OUTDIR}"/bayesprism_*.csv "${OUTDIR}"/bayesprism_*.rds; do
  if [[ -f "$f" ]]; then
    ext="${f##*.}"
    base="${f%.*}"
    # Handle .csv.gz double extension
    if [[ "$f" == *.csv.gz ]]; then
      backup="${f%.csv.gz}_star_backup.csv.gz"
    else
      backup="${base}_star_backup.${ext}"
    fi
    if [[ ! -f "$backup" ]]; then
      cp "$f" "$backup"
      echo "Backed up: $(basename $f) -> $(basename $backup)"
    fi
  fi
done

Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/118b_bayesprism_hires.R

echo ""
echo "Finished: $(date)"
echo "Outputs:"
ls -lh "${OUTDIR}"/bayesprism_* 2>/dev/null | grep -v star_backup
