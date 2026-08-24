#!/usr/bin/env bash
#SBATCH --job-name=hac_ml_pathway
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_ML_OUT_ROOT:?HAC_ML_OUT_ROOT is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

mapfile -t scripts < <(find "${SCRIPT_DIR}" -maxdepth 1 -type f \( -name '20_pathway_*.R' -o -name '21_tf_*.R' \) | sort)
[[ "${#scripts[@]}" -gt 0 ]] || { echo "No pathway/TF scripts found" >&2; exit 1; }
for script in "${scripts[@]}"; do "${R_BIN}" --vanilla "${script}"; done
echo "HAC_MOLECULAR_PATHWAY_TF_COMPLETE ${HAC_ML_OUT_ROOT}"
