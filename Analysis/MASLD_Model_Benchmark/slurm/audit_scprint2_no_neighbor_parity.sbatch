#!/usr/bin/env bash
#SBATCH --job-name=model-cpu-audit-302
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/scprint2-native-torch2.8.0-v2-21073597/source/.venv
SOURCE=${ROOT}/executions/transcriptformer-activation-preflight-21070971/fixture/unlabeled_counts.h5ad
FIXTURE=${ROOT}/executions/scprint2-no-neighbor-fixture-21074469
CELL_EMBEDDER=${ROOT}/executions/environments/scprint2-native-torch2.8.0-v2-21073597/source/scprint2/tasks/cell_emb.py
TARGET=${ROOT}/executions/scprint2-no-neighbor-parity-audit-${SLURM_JOB_ID}

if [[ -e "${TARGET}" ]]; then
  printf 'Refusing to overwrite scPRINT-2 parity audit.\n' >&2
  exit 1
fi
cd "${ROOT}"
export PYTHONPATH=${ROOT}:${ROOT}/src
"${RUNTIME}/bin/python" -m unittest -v \
  tests.unit.test_audit_scprint2_no_neighbor_parity
"${RUNTIME}/bin/python" scripts/audit_scprint2_no_neighbor_parity.py \
  --source "${SOURCE}" \
  --fixture-root "${FIXTURE}" \
  --cell-embedder-source "${CELL_EMBEDDER}" \
  --output "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
