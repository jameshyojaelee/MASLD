#!/usr/bin/env bash
#SBATCH --job-name=model-data-403
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
CONTRACT="${ROOT}/config/evaluation/gse244832_label_free_atac_exchange.toml"
: "${VALIDATION_JOB_ID:?VALIDATION_JOB_ID must identify the afterok validation job}"
VALIDATION="${ROOT}/executions/model-check-333-${VALIDATION_JOB_ID}"
TARGET="${ROOT}/executions/gse244832-label-free-atac-exchange-axis-20260825"
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
  "402ade8ffd26c450f4a763591870868620692f96f97c92310dac9174570d02b5 config/evaluation/gse244832_label_free_atac_exchange.toml" \
  "748fc627208e6cdc9d2a3ba79b3330ada1407e4b67a074fbe8fe4bbce24f3fe4 scripts/build_gse244832_label_free_atac_exchange.py" \
  "0a9cef0db265043f4d0fad2282693eb4144214e8266f87ed78819a9143b47359 tests/unit/test_gse244832_label_free_atac_exchange.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${VALIDATION}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(root)
metadata = manifest.get("metadata", {})
if (
    metadata.get("artifact_class") != "gse244832_label_free_atac_exchange_validation"
    or metadata.get("condition_values_read") is not False
    or metadata.get("cross_assay_join_inferred") is not False
    or metadata.get("biological_outcome_arrays_read") is not False
    or metadata.get("champion_claim_allowed") is not False
    or metadata.get("status") != "passed"
):
    raise SystemExit("validation authority differs")
PY

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" \
  scripts/build_gse244832_label_free_atac_exchange.py \
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
    != "gse244832_label_free_atac_exchange_axis"
    or summary.get("donors") != 18
    or summary.get("membership_cells") != 88814
    or summary.get("donor_lineage_rows") != 90
    or summary.get("windows") != 32000
    or summary.get("rotations") != 2
    or summary.get("condition_values_read") is not False
    or summary.get("rna_assay_opened") is not False
    or summary.get("cross_assay_join_inferred") is not False
    or summary.get("biological_outcome_arrays_read") is not False
    or summary.get("metrics_calculated") is not False
    or summary.get("gse244832_champion_claim_allowed") is not False
):
    raise SystemExit("label-free exchange output differs")
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE244832 label-free ATAC exchange axis: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
