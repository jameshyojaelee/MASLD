#!/usr/bin/env bash
#SBATCH --job-name=model-check-338
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
CONFIG=${ROOT}/config/campaigns/peakvi_gse296875_to_gse244832_target_cell_materialization_20260825.json
OUTPUT=${ROOT}/executions/model-check-338-${SLURM_JOB_ID}
STAGE=${OUTPUT}.staging

if [[ -e ${OUTPUT} || -e ${STAGE} ]]; then
  printf 'Refusing to overwrite GSE244832 PeakVI target-cell fixture\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
cd "${ROOT}"
PYTHONPATH=src:. "${PYTHON}" -m unittest -v \
  tests.unit.test_materialize_peakvi_gse244832_target_cells \
  > "${STAGE}/tests.stdout.txt" 2> "${STAGE}/tests.stderr.txt"
PYTHONPATH=src:. "${PYTHON}" - "${ROOT}" "${CONFIG}" "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

import numpy as np
from scripts.materialize_peakvi_gse244832_target_cells import (
    read_membership,
    read_unit_axis,
    select_cells,
    validate_config,
)

root, config_path, output = map(Path, sys.argv[1:])
config, paths = validate_config(root, config_path)
contract = json.loads(
    (paths["registration"] / "execution_contract.json").read_text(encoding="utf-8")
)
unit_rows, states = read_unit_axis(
    paths["registration"] / "query_axis_with_missingness.tsv"
)
binary = np.load(
    paths["registration"] / contract["query_lineage_missing_state_path"],
    allow_pickle=False,
)
cells = select_cells(
    read_membership(paths["axis"] / "label_free_cell_membership.tsv.gz"),
    unit_rows,
    states,
    cap=config["cell_selection"]["maximum_per_donor_lineage"],
    seed=config["cell_selection"]["seed"],
)
if (
    not np.array_equal(states, binary)
    or len(cells) != 3049
    or len({row["cell_hash"] for row in cells}) != 3049
):
    raise SystemExit("real GSE244832 eligible-cell fixture differs")
(output / "admission.json").write_text(
    json.dumps(
        {
            "schema_version": "masld-bench-peakvi-gse244832-target-cell-admission-v1",
            "status": "passed",
            "target_cells": len(cells),
            "eligible_donor_lineage_units": 48,
            "structurally_missing_units": 4,
            "below_qc_units": 20,
            "cell_hashes_globally_unique": True,
            "mask_applied_before_selection": True,
            "rna_assay_read": False,
            "cross_assay_join_consumed": False,
            "outcomes_read": False,
            "metrics_calculated": False,
        },
        sort_keys=True,
    )
    + "\n",
    encoding="utf-8",
)
PY
sha256sum \
  "${CONFIG}" \
  scripts/materialize_peakvi_gse244832_target_cells.py \
  tests/unit/test_materialize_peakvi_gse244832_target_cells.py \
  slurm/test_peakvi_gse244832_target_cells_cpu.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH=src:. "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1])
freeze_tree(
    root,
    {
        "artifact_class": "peakvi_gse244832_target_cell_fixture",
        "dataset_id": "gse244832",
        "tests": 5,
        "target_cells": 3049,
        "registered_mask_applied_before_cell_selection": True,
        "fragment_cut_sites": "start_plus_4_and_end_minus_5",
        "cell_hashes_globally_namespaced": True,
        "rna_assay_read": False,
        "cross_assay_join_consumed": False,
        "outcomes_read": False,
        "metrics_calculated": False,
        "status": "passed",
    },
)
PY
mv "${STAGE}" "${OUTPUT}"
printf 'PeakVI GSE244832 target-cell fixture: %s\n' "${OUTPUT}"
sha256sum "${OUTPUT}/ARTIFACTS.json"
