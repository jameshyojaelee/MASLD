#!/usr/bin/env bash
#SBATCH --job-name=model-data-401
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
CONTRACT="${ROOT}/config/evaluation/gse281367_label_free_atac_exchange.toml"
VALIDATION="${ROOT}/executions/gse281367-label-free-atac-exchange-validation-21100454"
TARGET="${ROOT}/executions/gse281367-label-free-atac-exchange-axis-20260825"
STAGE="${TARGET}.staging-${SLURM_JOB_ID}"
FAILED="${TARGET}.failed-${SLURM_JOB_ID}"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite label-free exchange axis: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"

for pair in \
  "5329086965d81bfab7a616e1fa423754da614a11fbc95c452520aba3a0f42237 config/evaluation/gse281367_label_free_atac_exchange.toml" \
  "7c9e5dcc4de2d43e6785ef84599284be67ce7f52f30e94d492e264ed20c13b5f scripts/build_gse281367_label_free_atac_exchange.py" \
  "3705dfbf632b7585649c1b78c3d26cc8682d07156bcdd1534300ac55cf0f2dd1 tests/unit/test_gse281367_label_free_atac_exchange.py" \
  "30660405c45fa117206e23841ba2d1a7f7763f2a40087592b3922b50749faed5 executions/gse281367-label-free-atac-exchange-validation-21100454/ARTIFACTS.json"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" \
  scripts/build_gse281367_label_free_atac_exchange.py \
  --root "${ROOT}" \
  --contract "${CONTRACT}" \
  --output "${STAGE}" \
  --hash-workers 4

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(root)
summary = json.loads((root / "exchange_summary.json").read_text(encoding="utf-8"))
if (
    manifest.get("metadata", {}).get("artifact_class")
    != "gse281367_label_free_atac_exchange_axis"
    or manifest.get("metadata", {}).get("condition_values_read") is not False
    or manifest.get("metadata", {}).get("biological_outcome_arrays_read") is not False
    or summary.get("donors") != 12
    or summary.get("membership_cells") != 226224
    or summary.get("donor_lineage_rows") != 60
    or summary.get("windows") != 32000
    or summary.get("rotations") != 2
    or summary.get("condition_values_read") is not False
    or summary.get("biological_outcome_arrays_read") is not False
    or summary.get("metrics_calculated") is not False
    or summary.get("gse281367_champion_claim_allowed") is not False
):
    raise SystemExit("label-free exchange output differs")
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Label-free ATAC exchange axis: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
