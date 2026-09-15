#!/usr/bin/env bash
#SBATCH --job-name=model-data-881
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
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
RUNTIME_ARTIFACTS_SHA256=f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb
PYTHON=${RUNTIME}/env/bin/python
PARTICIPANTS=${ROOT}/executions/model-data-072-21080651-participants
PARTICIPANTS_ARTIFACTS_SHA256=6655f23b2092e163532113e74d75d3c2a93857ee870b4f0a94fdcb4fa7f9dbc0
SUMMARY=${ROOT}/executions/model-data-088-21112137-summarization/summary
MOLECULAR=${ROOT}/executions/model-data-880-21130257-gse135251-source/source/molecular
MOLECULAR_ARTIFACTS_SHA256=8bd9037c46660df56f48210498964349657dfc4ea60a88d5d3430be45e31f97a
TARGET=${ROOT}/executions/model-data-881-${SLURM_JOB_ID}-activation
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE83452 transfer-activation artifacts.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
install -d -m 0750 "${STAGE}/activation"
cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${ROOT}/src:${ROOT}"
printf '%s  %s\n' "${RUNTIME_ARTIFACTS_SHA256}" "${RUNTIME}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${PARTICIPANTS_ARTIFACTS_SHA256}" "${PARTICIPANTS}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${MOLECULAR_ARTIFACTS_SHA256}" "${MOLECULAR}/ARTIFACTS.json" \
  | sha256sum --check --strict
bash -n slurm/build_gse83452_baseline_nash_transfer_activation_cpu.sbatch
"${PYTHON}" - "${PARTICIPANTS}" "${SUMMARY}" "${MOLECULAR}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("input frozen-tree verification passed")
PY
"${PYTHON}" -m unittest -v \
  tests.unit.test_build_gse83452_participant_masks \
  tests.unit.test_build_gse83452_baseline_selection \
  tests.unit.test_build_gse83452_transfer_gene_axis \
  > "${STAGE}/unit_tests.txt" 2>&1
# Rebuild the label-blind mask from the CURRENT participant table.  The
# superseded model-data-067 table carries the retired gse83452_longitudinal_array
# family; the builder asserts the family, so a rebuild against the wrong table
# cannot pass silently.
"${PYTHON}" scripts/build_gse83452_participant_masks.py \
  --records "${PARTICIPANTS}/gse83452_records.tsv" \
  --summarized-columns "${SUMMARY}/summarized_columns.txt" \
  --output "${STAGE}/activation/participant_masks.json" \
  > "${STAGE}/mask_stdout.json" 2> "${STAGE}/mask_stderr.txt"
# The summarization run already embedded a mask built by the same script.  A
# byte-identical rebuild is an independent check that the frozen copy was not
# produced from the superseded participant table.
cmp "${STAGE}/activation/participant_masks.json" "${SUMMARY}/participant_masks.json"
"${PYTHON}" scripts/build_gse83452_baseline_selection.py \
  --masks "${STAGE}/activation/participant_masks.json" \
  --output "${STAGE}/activation/baseline_participant_masks.json" \
  --output-columns "${STAGE}/activation/baseline_columns.txt" \
  > "${STAGE}/baseline_selection_stdout.json" 2> "${STAGE}/baseline_selection_stderr.txt"
"${PYTHON}" scripts/build_gse83452_transfer_gene_axis.py \
  --source-feature-axis "${MOLECULAR}/rna_feature_axis.tsv" \
  --external-gene-matrix "${SUMMARY}/gse83452_gpl16686_gene_matrix.tsv" \
  --external-key-field ensembl_gene_id \
  --expected-source-genes 73720 \
  --expected-external-genes 21787 \
  --expected-external-columns 231 \
  --output-axis "${STAGE}/activation/common_stable_gene_axis.tsv" \
  --output-receipt "${STAGE}/activation/gene_axis_receipt.json" \
  > "${STAGE}/gene_axis_stdout.json" 2> "${STAGE}/gene_axis_stderr.txt"
"${PYTHON}" -m pip freeze --all > "${STAGE}/pip_freeze.txt"
sha256sum \
  scripts/build_gse83452_participant_masks.py \
  scripts/build_gse83452_baseline_selection.py \
  scripts/build_gse83452_transfer_gene_axis.py \
  tests/unit/test_build_gse83452_participant_masks.py \
  tests/unit/test_build_gse83452_baseline_selection.py \
  tests/unit/test_build_gse83452_transfer_gene_axis.py \
  slurm/build_gse83452_baseline_nash_transfer_activation_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  "${PARTICIPANTS}/ARTIFACTS.json" \
  "${SUMMARY}/ARTIFACTS.json" \
  "${MOLECULAR}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import csv
import json
import sys

from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
activation = root / "activation"
masks = json.loads((activation / "participant_masks.json").read_text(encoding="utf-8"))
baseline = json.loads(
    (activation / "baseline_participant_masks.json").read_text(encoding="utf-8")
)
axis_receipt = json.loads((activation / "gene_axis_receipt.json").read_text(encoding="utf-8"))
with (activation / "common_stable_gene_axis.tsv").open(encoding="utf-8", newline="") as handle:
    axis = list(csv.DictReader(handle, delimiter="\t"))
gene_ids = [row["stable_gene_id"] for row in axis]
columns = [
    line.strip()
    for line in (activation / "baseline_columns.txt").read_text(encoding="utf-8").splitlines()
    if line.strip()
]
if (
    masks.get("status") != "pass_label_blind_participant_masks"
    or masks.get("records") != 231
    or masks.get("participants") != 171
    or masks.get("labels_read") is not False
    or masks.get("outcome_columns_read") is not False
    or masks.get("withheld_fields_imputed") is not False
):
    raise SystemExit("GSE83452 participant-mask verification differs")
if (
    baseline.get("status") != "pass_label_blind_participant_masks"
    or baseline.get("records") != 152
    or baseline.get("participants") != 152
    or baseline.get("distinct_outer_groups") != 152
    or baseline.get("followup_records_excluded") != 79
    or baseline.get("arrays_per_participant_after_selection") != 1
    or baseline.get("labels_read") is not False
    or len(baseline.get("mask_records", [])) != 152
    or len(columns) != 152
    or columns != [record["sample_accession"] for record in baseline["mask_records"]]
):
    raise SystemExit("GSE83452 baseline-selection verification differs")
forbidden = {
    "fibrosis_stage_group",
    "primary_transfer_evaluable",
    "paired_expression_stress_evaluable",
    "paired_nash_transition_evaluable",
}
for record in baseline["mask_records"]:
    if forbidden & set(record):
        raise SystemExit("GSE83452 baseline mask leaked an outcome column")
    for field in ("nash_status", "age", "sex", "intervention"):
        if record[field] is not None or record[f"{field}_observed"] is not False:
            raise SystemExit("structurally withheld field was filled")
if (
    axis_receipt.get("status") != "pass_label_blind_shared_gene_axis"
    or axis_receipt.get("external_expression_values_read") is not False
    or axis_receipt.get("external_labels_read") is not False
    or axis_receipt.get("shared_genes") != len(gene_ids)
    or gene_ids != sorted(gene_ids)
    or len(set(gene_ids)) != len(gene_ids)
):
    raise SystemExit("GSE83452 shared gene-axis verification differs")
freeze_tree(
    activation,
    {
        "artifact_class": "gse83452_baseline_nash_transfer_activation",
        "records": 231,
        "baseline_participants": 152,
        "shared_genes": len(gene_ids),
        "external_labels_read": False,
        "external_expression_values_read": False,
        "outcome_columns_read": False,
        "status": "passed",
    },
)
verify_frozen_tree(activation)
freeze_tree(
    root,
    {
        "artifact_class": "gse83452_baseline_nash_transfer_activation_campaign",
        "external_labels_read": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
verify_frozen_tree(root)
print(json.dumps({"shared_genes": len(gene_ids), "baseline_participants": 152}, sort_keys=True))
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/activation/ARTIFACTS.json"
