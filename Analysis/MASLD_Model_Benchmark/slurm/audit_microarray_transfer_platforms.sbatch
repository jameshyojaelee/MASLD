#!/usr/bin/env bash
#SBATCH --job-name=model-data-067
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/model-data-062-21079818
SOURCE_ARTIFACTS_SHA256=4c1460c573285fa53fe6464f23498e9c27c610c4be14cf3df3cc2944282ae96d
PARTICIPANTS_V1=${ROOT}/executions/model-data-063-21079860
PARTICIPANTS_V1_ARTIFACTS_SHA256=b6846c72d85133990bfa3852d939a1383471b07b892534ef292aaea3af569ff9
RAW=${ROOT}/executions/model-data-065-21079931
RAW_ARTIFACTS_SHA256=8b30a2265c421f2917ed886d99cfaa08ee8d874abd40e68c8c72d48c726781d1
PARTICIPANTS=${ROOT}/executions/model-data-067-${SLURM_JOB_ID}-participants
PARTICIPANTS_STAGE=${PARTICIPANTS}.staging
PARTICIPANTS_FAILED=${PARTICIPANTS}.failed
PLATFORM=${ROOT}/executions/model-data-067-${SLURM_JOB_ID}-platform
PLATFORM_STAGE=${PLATFORM}.staging
PLATFORM_FAILED=${PLATFORM}.failed

preserve_failed_stages() {
  if [[ -d "${PARTICIPANTS_STAGE}" && ! -e "${PARTICIPANTS_FAILED}" ]]; then
    mv "${PARTICIPANTS_STAGE}" "${PARTICIPANTS_FAILED}"
  fi
  if [[ -d "${PLATFORM_STAGE}" && ! -e "${PLATFORM_FAILED}" ]]; then
    mv "${PLATFORM_STAGE}" "${PLATFORM_FAILED}"
  fi
}
trap preserve_failed_stages EXIT
for path in \
  "${PARTICIPANTS}" "${PARTICIPANTS_STAGE}" "${PARTICIPANTS_FAILED}" \
  "${PLATFORM}" "${PLATFORM_STAGE}" "${PLATFORM_FAILED}"; do
  if [[ -e "${path}" ]]; then
    printf 'Refusing to overwrite microarray platform audit artifact: %s\n' "${path}" >&2
    exit 1
  fi
done

cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
printf '%s  %s\n' "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PARTICIPANTS_V1_ARTIFACTS_SHA256}" "${PARTICIPANTS_V1}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${RAW_ARTIFACTS_SHA256}" "${RAW}/ARTIFACTS.json" \
  | sha256sum --check --strict
python - "${SOURCE}" "${PARTICIPANTS_V1}" "${RAW}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
python -m unittest -v \
  tests.unit.test_audit_microarray_transfer_sources \
  tests.unit.test_audit_microarray_transfer_participants \
  tests.unit.test_acquire_audit_microarray_transfer_raw \
  tests.unit.test_audit_microarray_transfer_platforms \
  tests.contract.test_microarray_transfer_contract

python scripts/audit_microarray_transfer_participants.py \
  --source "${SOURCE}" \
  --source-artifacts-sha256 "${SOURCE_ARTIFACTS_SHA256}" \
  --output "${PARTICIPANTS_STAGE}"
sha256sum \
  scripts/audit_microarray_transfer_participants.py \
  tests/unit/test_audit_microarray_transfer_participants.py \
  "${SOURCE}/ARTIFACTS.json" \
  > "${PARTICIPANTS_STAGE}/source.sha256"
python - "${PARTICIPANTS_STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "participant_audit.json").read_text(encoding="utf-8"))
if (
    audit.get("schema_version")
    != "masld-bench-microarray-transfer-participant-audit-v2"
    or audit["gse83452"].get("paired_expression_stress_evaluable") != 60
    or audit["gse83452"].get("paired_nash_transition_evaluable") != 54
    or audit.get("model_training_activated") is not False
):
    raise SystemExit("microarray participant v2 firewall differs")
freeze_tree(
    root,
    {
        "artifact_class": "microarray_transfer_participant_audit_v2",
        "cohort_family_ids": [
            "gse31803_gse49541_fibrosis_array",
            "antwerp_inserm_shared",
        ],
        "participants": {"GSE49541": 72, "GSE83452": 171},
        "records": {"GSE49541": 72, "GSE83452": 231},
        "paired_expression_stress_evaluable": 60,
        "paired_nash_transition_evaluable": 54,
        "participant_join_authoritative": True,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "status": "passed_participant_topology_source_phenotypes_v2",
    },
)
PY
mv "${PARTICIPANTS_STAGE}" "${PARTICIPANTS}"
PARTICIPANTS_ARTIFACTS_SHA256=$(sha256sum "${PARTICIPANTS}/ARTIFACTS.json" | cut -d' ' -f1)

python scripts/audit_microarray_transfer_platforms.py \
  --source "${SOURCE}" \
  --source-artifacts-sha256 "${SOURCE_ARTIFACTS_SHA256}" \
  --participants "${PARTICIPANTS}" \
  --participants-artifacts-sha256 "${PARTICIPANTS_ARTIFACTS_SHA256}" \
  --raw "${RAW}" \
  --raw-artifacts-sha256 "${RAW_ARTIFACTS_SHA256}" \
  --output "${PLATFORM_STAGE}"
sha256sum \
  scripts/acquire_audit_microarray_transfer_raw.py \
  scripts/audit_microarray_transfer_platforms.py \
  tests/unit/test_acquire_audit_microarray_transfer_raw.py \
  tests/unit/test_audit_microarray_transfer_platforms.py \
  tests/contract/test_microarray_transfer_contract.py \
  config/evaluation/microarray_transfer_preprocessing.toml \
  config/evaluation/gse49541_fibrosis_transfer_task.toml \
  config/evaluation/gse83452_baseline_nash_transfer_task.toml \
  config/evaluation/gse83452_paired_intervention_stress_task.toml \
  slurm/audit_microarray_transfer_platforms.sbatch \
  "${SOURCE}/ARTIFACTS.json" \
  "${PARTICIPANTS}/ARTIFACTS.json" \
  "${RAW}/ARTIFACTS.json" \
  > "${PLATFORM_STAGE}/source.sha256"
python - "${PLATFORM_STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads(
    (root / "platform_rights_overlap_audit.json").read_text(encoding="utf-8")
)
if (
    audit.get("official_platform_axis_complete") is not True
    or audit.get("raw_CEL_format_complete") is not True
    or audit.get("rights_audit_complete") is not True
    or audit.get("platform_to_entrez_mapping_complete") is not False
    or audit.get("normalization_run") is not False
    or audit.get("model_training_activated") is not False
):
    raise SystemExit("microarray platform audit firewall differs")
freeze_tree(
    root,
    {
        "artifact_class": "microarray_transfer_platform_rights_overlap_audit",
        "cohort_family_ids": [
            "gse31803_gse49541_fibrosis_array",
            "antwerp_inserm_shared",
        ],
        "platforms": {"GSE49541": "GPL570", "GSE83452": "GPL16686"},
        "CEL_records": {"GSE49541": 72, "GSE83452": 231},
        "CEL_format": "calvin_v1",
        "official_platform_axis_complete": True,
        "rights_audit_complete": True,
        "inductive_preprocessing_complete": False,
        "model_training_activated": False,
        "sealed_outcomes_read": False,
        "status": "passed_platform_axes_CEL_signatures_rights_and_overlap",
    },
)
PY
mv "${PLATFORM_STAGE}" "${PLATFORM}"
trap - EXIT
sha256sum "${PARTICIPANTS}/ARTIFACTS.json" "${PLATFORM}/ARTIFACTS.json"
