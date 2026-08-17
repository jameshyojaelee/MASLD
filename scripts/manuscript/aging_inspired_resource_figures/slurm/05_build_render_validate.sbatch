#!/bin/bash
#SBATCH --job-name=agingfig456
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/candidates/aging-inspired-resource-rebuild-candidate-2026-08-15-v8/logs/build_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/figures/candidates/aging-inspired-resource-rebuild-candidate-2026-08-15-v8/logs/build_%j.log

set -euo pipefail
umask 027

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export FIG_CAND_ROOT=${MASLD_PROJECT_ROOT}/figures/candidates/aging-inspired-resource-rebuild-candidate-2026-08-15-v8
export SCRIPT_ROOT=${MASLD_PROJECT_ROOT}/scripts/manuscript/aging_inspired_resource_figures
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PYTHONNOUSERSITE=1
python_bin=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python

echo "[figure] host=$(hostname) job=${SLURM_JOB_ID} start=$(date -u)"
"${python_bin}" "${SCRIPT_ROOT}/05_build_figure_sources.py"
micromamba run -n decontx Rscript "${SCRIPT_ROOT}/08_render_compact_panels.R"
"${python_bin}" -m pip freeze > "${FIG_CAND_ROOT}/pip_freeze.txt"
"${python_bin}" "${SCRIPT_ROOT}/07_validate_figure_candidate.py"
echo "[figure] done=$(date -u)"
