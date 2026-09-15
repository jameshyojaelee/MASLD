#!/bin/bash
#SBATCH --job-name=cam_final
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
set -euo pipefail
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export OMP_NUM_THREADS=1
cd "$MASLD_PROJECT_ROOT"
R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
S=scripts/analysis/cross_assay_modules
"$R" --vanilla "$S/15_integrate_states.R"
"$R" --vanilla "$S/16_render_panel.R"
"$PY" "$S/17_validate_candidate.py" --out-root "$CAM_OUT_ROOT"
