#!/usr/bin/env bash
#SBATCH --job-name=model-check-261
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME="${ROOT}/executions/environments/corgi-nvidia-pytorch-24.01-py3-recovery-21064690"
IMAGE="${RUNTIME}/nvidia-pytorch-24.01-py3.sif"
TARGET="${ROOT}/executions/model-check-261-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi FiLM-head adapter test: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
for pair in \
  "01de778c17c67d92de9a2e242350e17443dbfe76b7c5a5f8bc447eae5408cd3c ${RUNTIME}/ARTIFACTS.json" \
  "eace72d1b831a4834d692134749bb82093c29930fb3a0258ae667b8acb840959 ${IMAGE}" \
  "6dd6242696f458efbec72bf73a36938802487a6d0e8f4e1ae2a401a3872129ca ${ROOT}/src/masld_bench/adapters/corgi_film_head.py" \
  "cb0dbb46d2eb96b7ea83f054dea00c3536e8a824f3d4ae064b243167f4f47a57 ${ROOT}/tests/unit/test_corgi_film_head_adapter.py"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done
module purge
module load singularity/4.2.2
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" \
  --pwd "${ROOT}" \
  --env "PYTHONPATH=${ROOT}/src:${ROOT}" \
  --env PYTHONNOUSERSITE=1 --env PYTHONDONTWRITEBYTECODE=1 \
  --env PYTHONHASHSEED=20260825 \
  "${IMAGE}" python -m unittest -v tests.unit.test_corgi_film_head_adapter \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
singularity exec --cleanenv --containall --net --network none \
  --bind "${ROOT}:${ROOT}:ro" --pwd "${ROOT}" \
  "${IMAGE}" python -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
sha256sum \
  src/masld_bench/adapters/corgi_film_head.py \
  tests/unit/test_corgi_film_head_adapter.py \
  slurm/test_corgi_regular_film_head_adapter_cpu.sbatch \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
stderr = (root / "validation/tests.stderr.txt").read_text(encoding="utf-8")
if "Ran 5 tests" not in stderr or not stderr.rstrip().endswith("OK"):
    raise SystemExit("Corgi FiLM-head adapter unit-test receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "corgi_regular_film_plus_head_adapter_unit_validation",
        "model_id": "corgi_regular",
        "task_id": "rna_conditioned_atac",
        "adaptation_rung": "film_plus_head",
        "tests": 5,
        "development_outcomes_read": False,
        "histology_or_disease_labels_read": False,
        "test_or_sealed_features_or_outcomes_read": False,
        "checkpoint_loaded": False,
        "training_performed": False,
        "prediction_performed": False,
        "global_census_modified": False,
        "global_promotion_gate_modified": False,
        "status": "passed",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi FiLM-plus-head adapter unit validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
