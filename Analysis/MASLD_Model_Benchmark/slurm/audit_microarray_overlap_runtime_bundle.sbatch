#!/usr/bin/env bash
#SBATCH --job-name=model-data-068
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
NUMPY_PYTHON=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ROOT=${PROJECT}/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/model-data-062-21079818
SOURCE_ARTIFACTS_SHA256=4c1460c573285fa53fe6464f23498e9c27c610c4be14cf3df3cc2944282ae96d
RAW=${ROOT}/executions/model-data-065-21079931
RAW_ARTIFACTS_SHA256=8b30a2265c421f2917ed886d99cfaa08ee8d874abd40e68c8c72d48c726781d1
PLATFORM_V1=${ROOT}/executions/model-data-067-21080193-platform
PLATFORM_V1_ARTIFACTS_SHA256=0ce6af4196d2fb0905d6ca3d07fbbeb3bf1ce5b37776a63a726cc1d21794f4bf
GSE106737_MATRIX=${PROJECT}/Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/sources/GSE106737/GSE106737_series_matrix.txt.gz
GSE106737_MATRIX_SHA256=d6c54bdb441f53f9fccc6d9592413a5c51a5cdce34900bb0b37ea9099e52cb7d
GSE83452_MATRIX=${PROJECT}/Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/sources/GSE83452/GSE83452_series_matrix.txt.gz
GSE83452_MATRIX_SHA256=f704bf623e4714e0c38317cc9bdf37e064ddae41185d1485b079e6fdc224603d
GSE106737_MANIFEST=${PROJECT}/Analysis/Multimodal_Program_Projection/candidates/public-functional-map-2026-08-09/sample_manifest.tsv
GSE106737_MANIFEST_SHA256=d71d905da5ed1eecae39e5d786f1065b02044c83c51404567b4d555dc06bc824
PRIOR_FINGERPRINT=${PROJECT}/Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09/source_gates/sample_fingerprint_matches.tsv
PRIOR_FINGERPRINT_SHA256=d13e209f4e6b93bc1de7e97340b6b05fcbaa5af65b69a61d2e8b23cf10854829
TIMEPOINT_DIAGNOSTIC=${ROOT}/executions/model-data-071-21080632-timepoints
TIMEPOINT_DIAGNOSTIC_ARTIFACTS_SHA256=53b625fd4cb514929b789d0e88be990600973d92924c15907da2cf7d1062b4bf

RUN_STEM=${SLURM_JOB_NAME}-${SLURM_JOB_ID}
PARTICIPANTS=${ROOT}/executions/${RUN_STEM}-participants
OVERLAP=${ROOT}/executions/${RUN_STEM}-overlap
PLATFORM=${ROOT}/executions/${RUN_STEM}-platform
PREPARATION=${ROOT}/executions/${RUN_STEM}-preparation
INDUCTIVE=${ROOT}/executions/${RUN_STEM}-inductive
RUNTIME=${ROOT}/executions/${RUN_STEM}-runtime
PARTICIPANTS_STAGE=${PARTICIPANTS}.staging
OVERLAP_STAGE=${OVERLAP}.staging
PLATFORM_STAGE=${PLATFORM}.staging
PREPARATION_STAGE=${PREPARATION}.staging
INDUCTIVE_STAGE=${INDUCTIVE}.staging
RUNTIME_STAGE=${RUNTIME}.staging
RECEIPTS_STAGE=${RUNTIME}.receipts.staging

preserve_failed_stages() {
  for stage in \
    "${PARTICIPANTS_STAGE}" "${OVERLAP_STAGE}" "${PLATFORM_STAGE}" \
    "${PREPARATION_STAGE}" "${INDUCTIVE_STAGE}" "${RUNTIME_STAGE}" \
    "${RECEIPTS_STAGE}"; do
    if [[ -e "${stage}" && ! -e "${stage}.failed" ]]; then
      mv "${stage}" "${stage}.failed"
    fi
  done
}
trap preserve_failed_stages EXIT
for path in \
  "${PARTICIPANTS}" "${OVERLAP}" "${PLATFORM}" "${PREPARATION}" \
  "${INDUCTIVE}" "${RUNTIME}" "${PARTICIPANTS_STAGE}" "${OVERLAP_STAGE}" \
  "${PLATFORM_STAGE}" "${PREPARATION_STAGE}" "${INDUCTIVE_STAGE}" \
  "${RUNTIME_STAGE}" "${RECEIPTS_STAGE}"; do
  if [[ -e "${path}" || -e "${path}.failed" ]]; then
    printf 'Refusing to overwrite microarray audit artifact: %s\n' "${path}" >&2
    exit 1
  fi
done

cd "${ROOT}"
module purge
module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
printf '%s  %s\n' \
  "${SOURCE_ARTIFACTS_SHA256}" "${SOURCE}/ARTIFACTS.json" \
  "${RAW_ARTIFACTS_SHA256}" "${RAW}/ARTIFACTS.json" \
  "${PLATFORM_V1_ARTIFACTS_SHA256}" "${PLATFORM_V1}/ARTIFACTS.json" \
  "${GSE106737_MATRIX_SHA256}" "${GSE106737_MATRIX}" \
  "${GSE83452_MATRIX_SHA256}" "${GSE83452_MATRIX}" \
  "${GSE106737_MANIFEST_SHA256}" "${GSE106737_MANIFEST}" \
  "${PRIOR_FINGERPRINT_SHA256}" "${PRIOR_FINGERPRINT}" \
  "${TIMEPOINT_DIAGNOSTIC_ARTIFACTS_SHA256}" "${TIMEPOINT_DIAGNOSTIC}/ARTIFACTS.json" \
  | sha256sum --check --strict
python - "${SOURCE}" "${RAW}" "${PLATFORM_V1}" "${TIMEPOINT_DIAGNOSTIC}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
python -m py_compile \
  scripts/audit_microarray_transfer_participants.py \
  scripts/audit_microarray_transfer_platforms.py \
  scripts/audit_microarray_cross_accession_overlap.py \
  scripts/prepare_microarray_runtime_audit.py \
  scripts/verify_microarray_inductive_fixture.py \
  scripts/finalize_microarray_runtime_audit.py
python -m unittest -v \
  tests.unit.test_audit_microarray_transfer_sources \
  tests.unit.test_audit_microarray_transfer_participants \
  tests.unit.test_acquire_audit_microarray_transfer_raw \
  tests.unit.test_audit_microarray_transfer_platforms \
  tests.contract.test_microarray_transfer_contract
"${NUMPY_PYTHON}" -m unittest -v \
  tests.unit.test_audit_microarray_cross_accession_overlap \
  tests.unit.test_verify_microarray_inductive_fixture

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
    audit["gse83452"].get("cohort_family_id") != "antwerp_inserm_shared"
    or audit["gse83452"].get("paired_expression_stress_evaluable") != 60
    or audit["gse83452"].get("paired_nash_transition_evaluable") != 54
    or audit.get("model_training_activated") is not False
):
    raise SystemExit("microarray participant v3 firewall differs")
freeze_tree(root, {
    "artifact_class": "microarray_transfer_participant_audit_v3",
    "cohort_family_ids": ["gse31803_gse49541_fibrosis_array", "antwerp_inserm_shared"],
    "participants": {"GSE49541": 72, "GSE83452": 171},
    "records": {"GSE49541": 72, "GSE83452": 231},
    "cross_accession_aliases_frozen_separately": True,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed_participant_topology_source_phenotypes_v3",
})
PY
mv "${PARTICIPANTS_STAGE}" "${PARTICIPANTS}"
PARTICIPANTS_ARTIFACTS_SHA256=$(sha256sum "${PARTICIPANTS}/ARTIFACTS.json" | cut -d' ' -f1)

"${NUMPY_PYTHON}" scripts/audit_microarray_cross_accession_overlap.py \
  --gse106737-series-matrix "${GSE106737_MATRIX}" \
  --gse83452-series-matrix "${GSE83452_MATRIX}" \
  --gse106737-manifest "${GSE106737_MANIFEST}" \
  --gse83452-records "${PARTICIPANTS}/gse83452_records.tsv" \
  --timepoint-diagnostic "${TIMEPOINT_DIAGNOSTIC}" \
  --timepoint-diagnostic-artifacts-sha256 "${TIMEPOINT_DIAGNOSTIC_ARTIFACTS_SHA256}" \
  --output "${OVERLAP_STAGE}"
python - "${OVERLAP_STAGE}" "${PRIOR_FINGERPRINT}" <<'PY'
from pathlib import Path
import csv
import sys
new_root, prior_path = map(Path, sys.argv[1:])
with (new_root / "cross_accession_sample_aliases.tsv").open(newline="", encoding="utf-8") as handle:
    new = {(row["gse106737_sample_accession"], row["gse83452_sample_accession"]) for row in csv.DictReader(handle, delimiter="\t")}
with prior_path.open(newline="", encoding="utf-8") as handle:
    prior = {(row["gse106737_sample"], row["gse83452_best_sample"]) for row in csv.DictReader(handle, delimiter="\t") if row["duplicate_threshold_pass"] == "TRUE"}
if len(new) != 78 or new != prior:
    raise SystemExit("independent overlap re-audit differs from prior frozen evidence")
PY
sha256sum \
  scripts/audit_microarray_cross_accession_overlap.py \
  tests/unit/test_audit_microarray_cross_accession_overlap.py \
  "${PARTICIPANTS}/ARTIFACTS.json" \
  "${GSE106737_MATRIX}" "${GSE83452_MATRIX}" "${GSE106737_MANIFEST}" \
  "${PRIOR_FINGERPRINT}" "${TIMEPOINT_DIAGNOSTIC}/ARTIFACTS.json" \
  > "${OVERLAP_STAGE}/source.sha256"
python - "${OVERLAP_STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "cross_accession_overlap_audit.json").read_text(encoding="utf-8"))
if (
    audit.get("cohort_family_id") != "antwerp_inserm_shared"
    or audit.get("reciprocal_array_aliases") != 78
    or audit.get("cross_accession_participant_aliases") != 41
    or audit.get("labels_read_for_overlap_audit") is not False
    or audit.get("expression_values_exported") is not False
    or audit.get("source_timepoints_retained_verbatim") is not True
    or audit.get("timepoints_coerced_or_overwritten") is not False
):
    raise SystemExit("microarray cross-accession overlap firewall differs")
freeze_tree(root, {
    "artifact_class": "microarray_cross_accession_overlap_audit",
    "cohort_family_id": "antwerp_inserm_shared",
    "accession_aliases": ["GSE106737", "GSE83452"],
    "reciprocal_array_aliases": 78,
    "cross_accession_participant_aliases": 41,
    "family_union_participant_range": [171, 200],
    "series_matrix_role": "contamination_audit_only_prohibited_model_input",
    "timepoint_diagnostic_artifacts_sha256": "53b625fd4cb514929b789d0e88be990600973d92924c15907da2cf7d1062b4bf",
    "source_timepoints_retained_verbatim": True,
    "timepoints_coerced_or_overwritten": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed_confirmed_cross_accession_family_overlap",
})
PY
mv "${OVERLAP_STAGE}" "${OVERLAP}"
OVERLAP_ARTIFACTS_SHA256=$(sha256sum "${OVERLAP}/ARTIFACTS.json" | cut -d' ' -f1)

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
  tests/unit/test_audit_microarray_transfer_platforms.py \
  "${SOURCE}/ARTIFACTS.json" "${PARTICIPANTS}/ARTIFACTS.json" \
  "${RAW}/ARTIFACTS.json" "${OVERLAP}/ARTIFACTS.json" \
  > "${PLATFORM_STAGE}/source.sha256"
python - "${PLATFORM_STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "platform_rights_overlap_audit.json").read_text(encoding="utf-8"))
if (
    audit["cohort_overlap"]["GSE83452"].get("cohort_family_id") != "antwerp_inserm_shared"
    or audit.get("platform_to_entrez_mapping_complete") is not False
    or audit.get("normalization_run") is not False
):
    raise SystemExit("microarray platform v2 firewall differs")
freeze_tree(root, {
    "artifact_class": "microarray_transfer_platform_rights_overlap_audit_v2",
    "cohort_family_ids": ["gse31803_gse49541_fibrosis_array", "antwerp_inserm_shared"],
    "cross_accession_overlap_artifacts_sha256": "provided_in_source_sha256",
    "platforms": {"GSE49541": "GPL570", "GSE83452": "GPL16686"},
    "official_platform_axis_complete": True,
    "rights_audit_complete": True,
    "inductive_preprocessing_complete": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed_platform_axes_CEL_signatures_rights_and_family_aliases_v2",
})
PY
mv "${PLATFORM_STAGE}" "${PLATFORM}"
PLATFORM_ARTIFACTS_SHA256=$(sha256sum "${PLATFORM}/ARTIFACTS.json" | cut -d' ' -f1)

python scripts/prepare_microarray_runtime_audit.py \
  --raw "${RAW}" \
  --raw-artifacts-sha256 "${RAW_ARTIFACTS_SHA256}" \
  --platform "${PLATFORM}" \
  --platform-artifacts-sha256 "${PLATFORM_ARTIFACTS_SHA256}" \
  --output "${PREPARATION_STAGE}"
sha256sum \
  scripts/prepare_microarray_runtime_audit.py \
  "${RAW}/ARTIFACTS.json" "${PLATFORM}/ARTIFACTS.json" \
  > "${PREPARATION_STAGE}/source.sha256"
python - "${PREPARATION_STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "runtime_preparation_audit.json").read_text(encoding="utf-8"))
if (
    audit.get("gpl16686_to_gencode_v49_crosswalk_complete") is not False
    or audit.get("gene_mapping_inferred_from_GB_ACC_or_coordinates") is not False
    or audit.get("normalization_run") is not False
):
    raise SystemExit("microarray runtime preparation firewall differs")
freeze_tree(root, {
    "artifact_class": "microarray_reader_fixture_and_annotation_disposition",
    "GPL16686_mapping_availability": audit["gpl16686_mapping_availability"],
    "biological_CEL_fixtures": 2,
    "normalization_run": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed_reader_fixtures_and_annotation_disposition",
})
PY
mv "${PREPARATION_STAGE}" "${PREPARATION}"

"${NUMPY_PYTHON}" scripts/verify_microarray_inductive_fixture.py --output "${INDUCTIVE_STAGE}"
sha256sum \
  scripts/verify_microarray_inductive_fixture.py \
  tests/unit/test_verify_microarray_inductive_fixture.py \
  config/evaluation/microarray_transfer_preprocessing.toml \
  > "${INDUCTIVE_STAGE}/source.sha256"
python - "${INDUCTIVE_STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "inductive_fixture_audit.json").read_text(encoding="utf-8"))
if (
    audit.get("tie_transform_invariant_to_feature_permutation") is not True
    or audit.get("tie_rank_invariant_to_feature_permutation") is not True
    or audit.get("biological_CEL_summarization_activated") is not False
):
    raise SystemExit("microarray inductive fixture firewall differs")
freeze_tree(root, {
    "artifact_class": "microarray_inductive_algorithmic_fixture",
    "fixture_role": "algorithmic_only_not_biological_preprocessing",
    "tie_safe": True,
    "held_array_isolation": True,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed_synthetic_held_array_isolation",
})
PY
mv "${INDUCTIVE_STAGE}" "${INDUCTIVE}"

mkdir "${RECEIPTS_STAGE}"
GPL570_FIXTURE=$(find "${PREPARATION}/fixtures" -maxdepth 1 -type f -name 'GSE49541_*.CEL' -print -quit)
GPL16686_FIXTURE=$(find "${PREPARATION}/fixtures" -maxdepth 1 -type f -name 'GSE83452_*.CEL' -print -quit)
module purge
module load R/4.4.3
Rscript scripts/audit_microarray_r_runtime.R \
  module_R_4.4.3 "${GPL570_FIXTURE}" "${GPL16686_FIXTURE}" \
  "${RECEIPTS_STAGE}/module_R_4.4.3.tsv"
module purge
module load R/4.4.1
Rscript scripts/audit_microarray_r_runtime.R \
  module_R_4.4.1 "${GPL570_FIXTURE}" "${GPL16686_FIXTURE}" \
  "${RECEIPTS_STAGE}/module_R_4.4.1.tsv"
module purge
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript \
  scripts/audit_microarray_r_runtime.R \
  micromamba_rnaseq "${GPL570_FIXTURE}" "${GPL16686_FIXTURE}" \
  "${RECEIPTS_STAGE}/micromamba_rnaseq.tsv"

module load python/3.11.5-GCCcore-13.2.0
export PYTHONPATH="${ROOT}:${ROOT}/src"
python scripts/finalize_microarray_runtime_audit.py \
  --preparation "${PREPARATION}" \
  --inductive-fixture "${INDUCTIVE}" \
  --runtime-receipt "${RECEIPTS_STAGE}/module_R_4.4.3.tsv" \
  --runtime-receipt "${RECEIPTS_STAGE}/module_R_4.4.1.tsv" \
  --runtime-receipt "${RECEIPTS_STAGE}/micromamba_rnaseq.tsv" \
  --preprocessing-contract config/evaluation/microarray_transfer_preprocessing.toml \
  --label-contract config/evaluation/microarray_transfer_label_semantics.toml \
  --output "${RUNTIME_STAGE}"
mv "${RECEIPTS_STAGE}"/*.tsv "${RUNTIME_STAGE}/"
rmdir "${RECEIPTS_STAGE}"
sha256sum \
  scripts/audit_microarray_r_runtime.R \
  scripts/finalize_microarray_runtime_audit.py \
  config/evaluation/microarray_transfer_preprocessing.toml \
  config/evaluation/microarray_transfer_label_semantics.toml \
  "${PREPARATION}/ARTIFACTS.json" "${INDUCTIVE}/ARTIFACTS.json" \
  > "${RUNTIME_STAGE}/source.sha256"
python - "${RUNTIME_STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "microarray_runtime_audit.json").read_text(encoding="utf-8"))
if (
    audit.get("synthetic_training_only_quantile_fixture_passed") is not True
    or audit.get("synthetic_fixture_is_biological_expression_authority") is not False
    or audit.get("GPL16686_to_GENCODE_v49_crosswalk_complete") is not False
    or audit.get("dependency_install_performed") is not False
    or audit.get("model_training_activated") is not False
):
    raise SystemExit("microarray runtime audit firewall differs")
freeze_tree(root, {
    "artifact_class": "microarray_existing_runtime_and_inductive_firewall_audit",
    "R_runtimes_audited": 3,
    "dependency_install_performed": False,
    "biological_preprocessing_activated": False,
    "model_training_activated": False,
    "sealed_outcomes_read": False,
    "status": "passed_runtime_inventory_biological_activation_blocked",
})
PY
mv "${RUNTIME_STAGE}" "${RUNTIME}"

trap - EXIT
sha256sum \
  "${PARTICIPANTS}/ARTIFACTS.json" \
  "${OVERLAP}/ARTIFACTS.json" \
  "${PLATFORM}/ARTIFACTS.json" \
  "${PREPARATION}/ARTIFACTS.json" \
  "${INDUCTIVE}/ARTIFACTS.json" \
  "${RUNTIME}/ARTIFACTS.json"
