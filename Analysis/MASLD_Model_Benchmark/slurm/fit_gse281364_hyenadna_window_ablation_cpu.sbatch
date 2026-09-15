#!/usr/bin/env bash
#SBATCH --job-name=hyenadna-window-ablation-fit
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ABLATION="${ROOT}/executions/gse281364-hyenadna-window-ablation-20260905T192600Z"
PYTHON="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
OUTPUT="${ABLATION}/evaluation"

install -d -m 0750 "${ABLATION}/validation"
if [[ -e "${OUTPUT}" ]]; then
  printf 'Refusing to overwrite window-ablation evaluation: %s\n' "${OUTPUT}" >&2
  exit 1
fi

export PYTHONPATH="${ROOT}:${ROOT}/src"
cd "${ROOT}"

"${PYTHON}" "${ROOT}/scripts/fit_gse281364_hyenadna_window_ablation.py" \
  --root "${ROOT}" \
  --ablation-embeddings "executions/gse281364-hyenadna-window-ablation-20260905T192600Z/embeddings" \
  --output "${OUTPUT}" \
  > "${ABLATION}/validation/fit.stdout.json" \
  2> "${ABLATION}/validation/fit.stderr.txt"

"${PYTHON}" -m pip freeze > "${ABLATION}/validation/fit.pip-freeze.txt"
"${PYTHON}" --version > "${ABLATION}/validation/fit.python.version.txt" 2>&1
sha256sum \
  "${ROOT}/scripts/fit_gse281364_hyenadna_window_ablation.py" \
  "${ROOT}/scripts/fit_gse281364_dna_language_seeded_heads.py" \
  "${ROOT}/scripts/fit_gse281364_open_sequence_heads.py" \
  "${ROOT}/scripts/gse281364_dna_lm_native_contract.py" \
  "${ROOT}/scripts/evaluate_gse281364_dna_lm_common_lane.py" \
  "${ABLATION}/prespec/PRESPECIFICATION.md" \
  > "${ABLATION}/validation/fit.source.sha256"

printf 'Completed HyenaDNA window-ablation fit and evaluation: %s\n' "${OUTPUT}"
