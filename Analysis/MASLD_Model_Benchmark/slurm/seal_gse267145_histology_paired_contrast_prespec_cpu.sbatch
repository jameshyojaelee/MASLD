#!/usr/bin/env bash
#SBATCH --job-name=model-check-1101
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
PRESPEC=${ROOT}/config/evaluation/gse267145_histology_paired_contrast_prespec_v1.json
PRESPEC_SHA256=fd00b3d2315fa9f3e729fb81e58ee3866c1f581c67d1aaed6a79c3b5cc7e094a
SEALER=${ROOT}/scripts/seal_gse267145_histology_paired_contrast_prespec.py
TARGET=${ROOT}/executions/model-check-1101-${SLURM_JOB_ID}-paired-contrast-prespec
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite the sealed paired-contrast prespecification.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
cd "${ROOT}"
export PYTHONPATH="${ROOT}/src:${ROOT}"

for pair in \
  "${PRESPEC_SHA256} ${PRESPEC}" \
  "dad1fa5b68686eaea6eeaa9c9725a4b0aed23dd4e96fcedbe1857d02ea070a43 ${ROOT}/scripts/score_gse267145_histology_production_v4.py" \
  "1d753b8e23fda8f9582e4d8cc2ce9c90eeae65d370ce342bf9d32ffdd2ba5920 ${ROOT}/scripts/evaluate_gse267145_histology_predictions.py" \
  "53a85a843e75b1c0fa088af86cb2f2593976d28039807dfa17d1ec571603778c ${ROOT}/config/evaluation/gse267145_histology_state_task.toml" \
  "532f66eab44561b4b84124c745b24a839f7d267c970fe5fb1fa9206513be7971 ${ROOT}/config/evaluation/gse267145_histology_production_v4_scoring.json" \
  "c8fd739e7d9cbaca9bf64729a8748dd0add42a9c76a5f7d4088362caf5add494 ${ROOT}/config/evaluation/gse267145_histology_state_promotion_gate.json" \
  "a55aed3d9ac4cad2211ed3b319c45f911d05851b8559408ed99bf7fa59562393 ${ROOT}/executions/model-scoring-078-21092130/scores/ARTIFACTS.json" \
  "98a74b97f6a9712a109d2516c3803a4acf15030ed607f861e1ce8c22b4f9753a ${ROOT}/executions/model-data-061-21079623/activation/ARTIFACTS.json" \
  "9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3 ${ROOT}/executions/model-data-064-21079902/fixture/folds/ARTIFACTS.json" \
  "f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb ${RUNTIME}/ARTIFACTS.json"; do
  read -r expected path <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${path}" | sha256sum --check --strict
done

"${PYTHON}" "${SEALER}" \
  --root "${ROOT}" \
  --prespec "${PRESPEC}" \
  --prespec-sha256 "${PRESPEC_SHA256}" \
  --output "${STAGE}/prespec" \
  > "${STAGE}/validation/sealer_stdout.json" \
  2> "${STAGE}/validation/sealer_stderr.txt"

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  "${PRESPEC}" \
  "${SEALER}" \
  "${ROOT}/scripts/analyze_gse267145_histology_paired_contrasts.py" \
  slurm/seal_gse267145_histology_paired_contrast_prespec_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/validation/source.sha256"

"${PYTHON}" - "${STAGE}" "${PRESPEC_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
prespec_sha256 = sys.argv[2]
verify_frozen_tree(root / "prespec")
receipt = json.loads(
    (root / "prespec/prespec_seal_receipt.json").read_text(encoding="utf-8")
)
if (
    receipt.get("prespec_contract_sha256") != prespec_sha256
    or receipt.get("declared_contrasts") != 189
    or receipt.get("computable_contrasts") != 169
    or receipt.get("not_applicable_contrasts") != 20
    or receipt.get("outcomes_read") is not False
    or receipt.get("metrics_calculated") is not False
    or receipt.get("contrasts_calculated") is not False
    or receipt.get("p_values_calculated") is not False
    or receipt.get("bh_adjustment_calculated") is not False
    or receipt.get("model_ranking_performed") is not False
):
    raise SystemExit("sealed prespecification receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "gse267145_histology_paired_contrast_prespec_campaign",
        "prespec_contract_sha256": prespec_sha256,
        "declared_contrasts": 189,
        "computable_contrasts": 169,
        "not_applicable_contrasts": 20,
        "outcomes_read": False,
        "metrics_calculated": False,
        "contrasts_calculated": False,
        "status": "sealed_prespec_no_outcome_read",
    },
)
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/prespec/ARTIFACTS.json"
