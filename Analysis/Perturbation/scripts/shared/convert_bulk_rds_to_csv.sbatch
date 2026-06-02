#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:0
#SBATCH --time=2:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=2
#SBATCH --job-name=bulk-csv
#SBATCH --output=Analysis/Perturbation/results/model_availability/logs/bulk_csv_convert_%j.log

# One-shot conversion of 5-cohort bulk .rds → CSV for Python loaders.
# Uses rnaseq env (which has data.table + R Bioc stack).

# Note: NO set -u before micromamba activate (activate.d scripts use unbound vars)
set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

export PYTHONNOUSERSITE=1

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u  # safe now

echo "=== R version ==="
R --version | head -1
echo ""
echo "=== Running converter ==="
Rscript Analysis/Perturbation/scripts/shared/convert_bulk_rds_to_csv.R

echo ""
echo "=== Output ==="
ls -lh Analysis/Perturbation/data/pseudobulks/ 2>&1
echo ""
echo "[done] bulk-csv conversion complete"
