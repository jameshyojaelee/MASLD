#!/bin/bash
#SBATCH --job-name=masld_donor_curve_panel
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_donor_curve_panel_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/masld_donor_curve_panel_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python
OUT=$PROJ/GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model/donor_curve_panel_${SLURM_JOB_ID}
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_donor_curve_prepare.py" manifest --out "$OUT"
"$PY" "$PROJ/scripts/analysis/alphagenome_campaign/model_donor_curve_prepare.py" paired --out "$OUT"
