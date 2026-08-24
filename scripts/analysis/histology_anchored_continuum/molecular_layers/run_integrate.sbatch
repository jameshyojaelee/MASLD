#!/usr/bin/env bash
#SBATCH --job-name=hac_ml_integrate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/histology_anchored_continuum/molecular_layers/slurm-%x-%j.out

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="${PROJECT_ROOT}/scripts/analysis/histology_anchored_continuum/molecular_layers"
R_BIN="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
: "${HAC_ML_OUT_ROOT:?HAC_ML_OUT_ROOT is unset}"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

for script in "${SCRIPT_DIR}"/*.R "${SCRIPT_DIR}"/tests/*.R; do
  "${R_BIN}" --vanilla -e "parse(file='${script}')" >/dev/null
done
"${R_BIN}" --vanilla "${SCRIPT_DIR}/tests/run_tests.R" 2>&1 |
  tee "${HAC_ML_OUT_ROOT}/integration_contract_tests.${SLURM_JOB_ID}.log"
"${R_BIN}" --vanilla "${SCRIPT_DIR}/50_integrate_gates.R"
"${R_BIN}" --vanilla "${SCRIPT_DIR}/60_render_figures.R"
python3 "${SCRIPT_DIR}/70_validate_candidate.py" "${HAC_ML_OUT_ROOT}"
python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py" 2>&1 |
  tee "${HAC_ML_OUT_ROOT}/validation/resource_scope.log"
echo "HAC_MOLECULAR_CANDIDATE_COMPLETE ${HAC_ML_OUT_ROOT}"
