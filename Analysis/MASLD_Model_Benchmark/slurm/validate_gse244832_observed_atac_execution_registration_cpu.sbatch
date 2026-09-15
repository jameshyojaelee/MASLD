#!/bin/bash
#SBATCH --job-name=model-check-337
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out

set -euo pipefail

ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark"
CONFIG="${ROOT}/config/evaluation/gse244832_observed_atac_execution_registration_20260825.json"
PYTHON="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
TARGET="${ROOT}/executions/model-check-337-${SLURM_JOB_ID}"

cd "${ROOT}"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

sha256sum -c <<'EOF'
1e3199e3483fa865b039012f7df6e02f5424a5fb35226accd2aa3963264ca9a9  config/evaluation/gse244832_observed_atac_execution_registration_20260825.json
908c85acfddf55b8ec24d3f9297ae19d9cce218163a112351555918d64c3ca19  scripts/register_gse244832_observed_atac_execution.py
a556051a0b4b3e6561fbe52101e9d8fcb456fbe04e6c53df092e158cd6ff3f10  tests/unit/test_register_gse244832_observed_atac_execution.py
29f913b8736012462dddc92d72304354ee9c4c92b1ec0941a601ef097e09bbee  scripts/commit_gse244832_atac_transport_predictions.py
23a21e7903f64c4b85986a2df59122ec1ddc97d1737251dab036de0e1519de69  tests/unit/test_commit_gse244832_atac_transport_predictions.py
d9f90c35d60da54ebadf61925c5ae7e1d3ca907e15237dfc2fb353d0b80a8d98  scripts/score_gse244832_atac_transport.py
8b958372fbd63e061bd87b1e5e68d92028d834cce4ccfcbe861192b7e3611d53  tests/unit/test_score_gse244832_atac_transport.py
ff1a506d5cca81382911355978097ecd49fcb1b2b48826e8e760a5f2ad1ad35e  executions/gse244832-reference-guarded-atac-exchange-axis-20260825/ARTIFACTS.json
801411d062f9d9f199c0468b6c8c6256263af85ea2805d75e3d166299d190742  executions/gse244832-label-free-atac-query-valid-context-20260825/ARTIFACTS.json
970992d3dc84f2d2c6c2c5357149579fa09b16051ddafb7eed0c1f2949a9ee8a  executions/gse244832-label-free-atac-query-test-context-20260825/ARTIFACTS.json
EOF

PYTHONPATH="src:." "${PYTHON}" -m unittest \
  tests.unit.test_build_gse244832_reference_guarded_atac_exchange \
  tests.unit.test_materialize_gse244832_label_free_query_atac \
  tests.unit.test_register_gse244832_observed_atac_execution \
  tests.unit.test_commit_gse244832_atac_transport_predictions \
  tests.unit.test_score_gse244832_atac_transport

PYTHONPATH="src:." "${PYTHON}" scripts/register_gse244832_observed_atac_execution.py \
  --root "${ROOT}" \
  --config "${CONFIG}" \
  --output "${TARGET}"

PYTHONPATH="src:." "${PYTHON}" - "${TARGET}" <<'PY'
import json
from pathlib import Path
import sys

import numpy as np

from masld_bench.artifacts import verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(root)
contract = json.loads((root / "execution_contract.json").read_text(encoding="utf-8"))
receipt = json.loads((root / "validation_receipt.json").read_text(encoding="utf-8"))
states = np.load(root / "query_lineage_missing_state.uint8.npy", allow_pickle=False)
assert manifest["metadata"]["artifact_class"] == "gse244832_observed_atac_execution_registration"
assert manifest["metadata"]["allowed_model_ids"] == [
    "lsi", "observed_atac_glm", "observed_atac_only", "peakvi"
]
assert contract["prediction_commit_required_before_outcome_open"] is True
assert contract["biological_outcome_arrays_read"] is False
assert contract["query_values_reusable_only_through_this_registration"] is True
assert receipt["eligible_units"] == 48
assert receipt["structurally_missing_units"] == 4
assert receipt["below_qc_units"] == 20
assert states.shape == (18, 4)
assert set(np.unique(states)) == {0, 1, 3}
print(json.dumps(receipt, sort_keys=True))
PY

sha256sum "${TARGET}/ARTIFACTS.json"
