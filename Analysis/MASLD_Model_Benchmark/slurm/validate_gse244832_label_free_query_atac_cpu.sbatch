#!/usr/bin/env bash
#SBATCH --job-name=model-check-335
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PYTHON="${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python"
TARGET="${ROOT}/executions/model-check-335-${SLURM_JOB_ID}"
if [[ -e "${TARGET}" ]]; then
  printf 'Refusing to overwrite query validation: %s\n' "${TARGET}" >&2
  exit 1
fi
mkdir -p "${TARGET}"
cd "${ROOT}"

for pair in \
  "dfbf1a2db5b9d0eb4768081e61145360616e932720b036eed7fcfd09bebb4788 config/campaigns/gse244832_label_free_query_atac_materialization_20260825.json" \
  "1c6a99e7eb69b8c358da804d1b00730502b00c7a7b29fa9df2f790d62458749f scripts/materialize_gse244832_label_free_query_atac.py" \
  "9db6a001dde8170d7eea214b7dbbddf7c3b20b317ae65ca1c2b56242de8e6c81 tests/unit/test_materialize_gse244832_label_free_query_atac.py" \
  "ff1a506d5cca81382911355978097ecd49fcb1b2b48826e8e760a5f2ad1ad35e executions/gse244832-reference-guarded-atac-exchange-axis-20260825/ARTIFACTS.json"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" -m unittest -v \
  tests.unit.test_materialize_gse244832_label_free_query_atac
"${PYTHON}" -m pip freeze > "${TARGET}/pip-freeze.txt"

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" - "${ROOT}" "${TARGET}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, write_json_exclusive
from scripts.materialize_gse244832_label_free_query_atac import validate_config

root = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2]).resolve(strict=True)
config_path = root / "config/campaigns/gse244832_label_free_query_atac_materialization_20260825.json"
config = json.loads(config_path.read_text(encoding="utf-8"))
validate_config(root, config)
receipt = {
    "schema_version": "masld-bench-gse244832-query-atac-validation-v1",
    "dataset_id": "gse244832",
    "donors": 18,
    "fragment_cut_sites": "start_plus_4_and_end_minus_5",
    "sequence_execution_authorized": False,
    "observed_atac_query_materialization_authorized": True,
    "condition_values_read": False,
    "rna_assay_opened": False,
    "cross_assay_join_inferred": False,
    "evaluator_outcome_artifact_read": False,
    "metrics_calculated": False,
    "champion_claim_allowed": False,
    "status": "passed",
}
write_json_exclusive(target / "validation_receipt.json", receipt)
freeze_tree(
    target,
    {
        "artifact_class": "gse244832_label_free_query_atac_validation",
        "dataset_id": "gse244832",
        "sequence_execution_authorized": False,
        "condition_values_read": False,
        "rna_assay_opened": False,
        "evaluator_outcome_artifact_read": False,
        "metrics_calculated": False,
        "champion_claim_allowed": False,
        "status": "passed",
    },
)
print(json.dumps(receipt, sort_keys=True))
PY

printf 'GSE244832 query-ATAC validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
