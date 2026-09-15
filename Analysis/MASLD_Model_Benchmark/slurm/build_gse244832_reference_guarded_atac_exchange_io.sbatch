#!/usr/bin/env bash
#SBATCH --job-name=model-data-404
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
GUARD="${ROOT}/config/evaluation/gse244832_atac_reference_guard_20260825.json"
: "${VALIDATION_JOB_ID:?VALIDATION_JOB_ID must identify the afterok validation job}"
VALIDATION="${ROOT}/executions/model-check-334-${VALIDATION_JOB_ID}"
TARGET="${ROOT}/executions/gse244832-reference-guarded-atac-exchange-axis-20260825"
STAGE="${TARGET}.staging-${SLURM_JOB_ID}"
FAILED="${TARGET}.failed-${SLURM_JOB_ID}"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite reference-guarded axis: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"

for pair in \
  "402ade8ffd26c450f4a763591870868620692f96f97c92310dac9174570d02b5 config/evaluation/gse244832_label_free_atac_exchange.toml" \
  "da94da73720a1e799e49edcfc261f37ef53cf6bcbb6abd4c92773d3e84fef368 config/evaluation/gse244832_atac_reference_guard_20260825.json" \
  "881391442717d4eccc71a3aa0e1861fa8d1c027f19796a992da660eed17ff3c2 scripts/build_gse244832_reference_guarded_atac_exchange.py" \
  "5ce29f4faa7735c143c1b5f0244bdf5a01a7ccd986f9c26013964c61c6992fe3 tests/unit/test_build_gse244832_reference_guarded_atac_exchange.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${VALIDATION}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
metadata = verify_frozen_tree(root).get("metadata", {})
if (
    metadata.get("artifact_class") != "gse244832_reference_guarded_atac_exchange_validation"
    or metadata.get("sequence_execution_authorized") is not False
    or metadata.get("observed_atac_query_materialization_authorized") is not True
    or metadata.get("biological_outcome_arrays_read") is not False
    or metadata.get("status") != "passed"
):
    raise SystemExit("reference-guard validation authority differs")
PY

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" \
  scripts/build_gse244832_reference_guarded_atac_exchange.py \
  --root "${ROOT}" \
  --contract "${CONTRACT}" \
  --guard "${GUARD}" \
  --output "${STAGE}" \
  --hash-workers 4

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
metadata = verify_frozen_tree(root).get("metadata", {})
summary = json.loads((root / "exchange_summary.json").read_text(encoding="utf-8"))
if (
    metadata.get("artifact_class") != "gse244832_reference_guarded_atac_exchange_axis"
    or summary.get("donors") != 18
    or summary.get("membership_cells") != 88814
    or summary.get("donor_lineage_rows") != 90
    or summary.get("family_rows") != 25
    or summary.get("non_observed_atac_executable_rows") != 0
    or summary.get("sequence_or_coordinate_family_execution_authorized") is not False
    or summary.get("observed_atac_query_materialization_authorized") is not True
    or summary.get("condition_values_read") is not False
    or summary.get("rna_assay_opened") is not False
    or summary.get("biological_outcome_arrays_read") is not False
    or summary.get("metrics_calculated") is not False
):
    raise SystemExit("reference-guarded exchange output differs")
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE244832 reference-guarded ATAC axis: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
