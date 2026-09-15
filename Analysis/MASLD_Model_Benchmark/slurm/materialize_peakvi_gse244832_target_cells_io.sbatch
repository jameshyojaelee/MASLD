#!/usr/bin/env bash
#SBATCH --job-name=model-data-339
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CONFIG=${ROOT}/config/campaigns/peakvi_gse296875_to_gse244832_target_cell_materialization_20260825.json
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
VALID=${ROOT}/executions/peakvi-gse244832-valid-context-cells-${SLURM_JOB_ID}
TEST=${ROOT}/executions/peakvi-gse244832-test-context-cells-${SLURM_JOB_ID}
WORK=${ROOT}/executions/peakvi-gse244832-cell-materialization-work-${SLURM_JOB_ID}
FAILED=${ROOT}/executions/peakvi-gse244832-cell-materialization-failed-${SLURM_JOB_ID}

preserve_failure() {
  if [[ -d ${WORK} && ! -e ${FAILED} && ! -e ${WORK}/COMPLETE ]]; then
    mv "${WORK}" "${FAILED}"
  fi
}
trap preserve_failure EXIT
if [[ -e ${VALID} || -e ${TEST} || -e ${WORK} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite GSE244832 PeakVI target-cell materialization\n' >&2
  exit 1
fi

cd "${ROOT}"
PYTHONPATH=src:. "${PYTHON}" scripts/materialize_peakvi_gse244832_target_cells.py \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --target-valid-output "${VALID}" \
  --target-test-output "${TEST}" \
  --work "${WORK}" \
  --workers 6

PYTHONPATH=src:. "${PYTHON}" - "${WORK}" "${VALID}" "${TEST}" <<'PY'
from hashlib import sha256
from pathlib import Path
import csv
import json
import sys

import numpy as np
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

work, valid, test = map(Path, sys.argv[1:])
receipt = json.loads((work / "receipt.json").read_text(encoding="utf-8"))
if (
    receipt.get("status") != "passed"
    or receipt.get("eligible_donor_lineage_units") != 48
    or receipt.get("structurally_missing_units") != 4
    or receipt.get("below_qc_units") != 20
    or receipt.get("registered_mask_applied_before_cell_selection") is not True
    or receipt.get("ineligible_cells_entered_transform") is not False
    or receipt.get("fragment_cut_sites") != "start_plus_4_and_end_minus_5"
    or receipt.get("rna_assay_read") is not False
    or receipt.get("cross_assay_join_consumed") is not False
    or receipt.get("metrics_calculated") is not False
):
    raise SystemExit("GSE244832 target-cell materialization receipt differs")
axes = []
for role, child in (("valid", valid), ("test", test)):
    manifest = verify_frozen_tree(child)
    contract = json.loads((child / "contract.json").read_text(encoding="utf-8"))
    with (child / "cell_axis.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    states = np.load(child / "donor_lineage_missing_state.uint8.npy", allow_pickle=False)
    if (
        manifest.get("metadata", {}).get("artifact_class")
        != "peakvi_gse244832_target_cell_context"
        or contract.get("context_role") != role
        or contract.get("eligible_donor_lineage_units") != 48
        or contract.get("registered_mask_applied_before_cell_selection") is not True
        or contract.get("ineligible_cells_entered_transform") is not False
        or contract.get("fragment_cut_sites") != "start_plus_4_and_end_minus_5"
        or states.shape != (18, 4)
        or {value: int((states == value).sum()) for value in (0, 1, 3)}
        != {0: 48, 1: 4, 3: 20}
        or len({row["cell_hash"] for row in rows}) != len(rows)
    ):
        raise SystemExit(f"GSE244832 target-cell artifact differs: {role}")
    axes.append((child / "cell_axis.tsv").read_bytes())
if axes[0] != axes[1]:
    raise SystemExit("rotation cell axes are not byte-identical")
freeze_tree(
    work,
    {
        "artifact_class": "peakvi_gse244832_target_cell_materialization",
        "source_dataset_id": "gse296875",
        "target_dataset_id": "gse244832",
        "eligible_donor_lineage_units": 48,
        "structurally_missing_units": 4,
        "below_qc_units": 20,
        "fragment_cut_sites": "start_plus_4_and_end_minus_5",
        "cell_hashes_globally_unique": True,
        "rotations_share_identical_cell_axis": True,
        "rna_assay_read": False,
        "cross_assay_join_consumed": False,
        "evaluator_outcomes_accessed": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
PY
trap - EXIT
printf 'GSE244832 PeakVI target-cell artifacts:\n'
for path in "${VALID}" "${TEST}" "${WORK}"; do
  printf '%s\n' "${path}"
  sha256sum "${path}/ARTIFACTS.json"
done
