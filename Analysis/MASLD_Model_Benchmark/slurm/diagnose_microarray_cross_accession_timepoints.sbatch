#!/usr/bin/env bash
#SBATCH --job-name=model-data-071
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT}/Analysis/MASLD_Model_Benchmark
PARTICIPANTS=${ROOT}/executions/model-data-068-21080575-participants
PARTICIPANTS_ARTIFACTS_SHA256=6655f23b2092e163532113e74d75d3c2a93857ee870b4f0a94fdcb4fa7f9dbc0
FINGERPRINT=${PROJECT}/Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/source_gates/sample_fingerprint_matches.tsv
FINGERPRINT_SHA256=d13e209f4e6b93bc1de7e97340b6b05fcbaa5af65b69a61d2e8b23cf10854829
GSE106737_MANIFEST=${PROJECT}/Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/sample_manifest.tsv
GSE106737_MANIFEST_SHA256=d71d905da5ed1eecae39e5d786f1065b02044c83c51404567b4d555dc06bc824
OUTPUT=${ROOT}/executions/${SLURM_JOB_NAME}-${SLURM_JOB_ID}-timepoints
STAGE=${OUTPUT}.staging
FAILED=${OUTPUT}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
for path in "${OUTPUT}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite timepoint diagnostic: %s\n' "${path}" >&2
    exit 1
  fi
done

cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
printf '%s  %s\n' \
  "${PARTICIPANTS_ARTIFACTS_SHA256}" "${PARTICIPANTS}/ARTIFACTS.json" \
  "${FINGERPRINT_SHA256}" "${FINGERPRINT}" \
  "${GSE106737_MANIFEST_SHA256}" "${GSE106737_MANIFEST}" \
  | sha256sum --check --strict
python - "${PARTICIPANTS}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
verify_frozen_tree(Path(sys.argv[1]))
PY
python -m py_compile scripts/diagnose_microarray_cross_accession_timepoints.py
python -m unittest -v tests.unit.test_diagnose_microarray_cross_accession_timepoints
python scripts/diagnose_microarray_cross_accession_timepoints.py \
  --fingerprint "${FINGERPRINT}" \
  --fingerprint-sha256 "${FINGERPRINT_SHA256}" \
  --gse106737-manifest "${GSE106737_MANIFEST}" \
  --gse106737-manifest-sha256 "${GSE106737_MANIFEST_SHA256}" \
  --gse83452-records "${PARTICIPANTS}/gse83452_records.tsv" \
  --output "${STAGE}"
sha256sum \
  scripts/diagnose_microarray_cross_accession_timepoints.py \
  tests/unit/test_diagnose_microarray_cross_accession_timepoints.py \
  "${PARTICIPANTS}/ARTIFACTS.json" "${FINGERPRINT}" "${GSE106737_MANIFEST}" \
  > "${STAGE}/source.sha256"
python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "cross_accession_timepoint_diagnostic.json").read_text(encoding="utf-8"))
if (
    audit.get("array_aliases_audited") != 78
    or audit.get("participant_aliases_audited") != 41
    or audit.get("source_timepoints_retained_verbatim") is not True
    or audit.get("timepoints_coerced_or_overwritten") is not False
    or audit.get("labels_read") is not False
):
    raise SystemExit("cross-accession timepoint diagnostic firewall differs")
freeze_tree(root, {
    "artifact_class": "microarray_cross_accession_timepoint_diagnostic",
    "cohort_family_id": "antwerp_inserm_shared",
    "array_aliases_audited": 78,
    "participant_aliases_audited": 41,
    "source_timepoints_retained_verbatim": True,
    "timepoints_coerced_or_overwritten": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": audit["status"],
})
PY
mv "${STAGE}" "${OUTPUT}"
trap - EXIT
sha256sum "${OUTPUT}/ARTIFACTS.json"
