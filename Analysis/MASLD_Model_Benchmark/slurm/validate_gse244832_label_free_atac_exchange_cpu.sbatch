#!/usr/bin/env bash
#SBATCH --job-name=model-check-333
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
TARGET="${ROOT}/executions/model-check-333-${SLURM_JOB_ID}"
if [[ -e "${TARGET}" ]]; then
  printf 'Refusing to overwrite validation artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
mkdir -p "${TARGET}"
cd "${ROOT}"

for pair in \
  "402ade8ffd26c450f4a763591870868620692f96f97c92310dac9174570d02b5 config/evaluation/gse244832_label_free_atac_exchange.toml" \
  "748fc627208e6cdc9d2a3ba79b3330ada1407e4b67a074fbe8fe4bbce24f3fe4 scripts/build_gse244832_label_free_atac_exchange.py" \
  "0a9cef0db265043f4d0fad2282693eb4144214e8266f87ed78819a9143b47359 tests/unit/test_gse244832_label_free_atac_exchange.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" -m unittest -v \
  tests.unit.test_gse244832_label_free_atac_exchange \
  tests.unit.test_gse281367_label_free_atac_exchange

"${PYTHON}" -m pip freeze > "${TARGET}/pip-freeze.txt"
PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" - "${ROOT}" "${TARGET}" <<'PY'
from __future__ import annotations

import json
from pathlib import Path
import sys
import tomllib

from masld_bench.artifacts import freeze_tree, write_json_exclusive
from scripts.build_gse244832_label_free_atac_exchange import validate_contract

root = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2]).resolve(strict=True)
contract_path = root / "config/evaluation/gse244832_label_free_atac_exchange.toml"
with contract_path.open("rb") as handle:
    contract = tomllib.load(handle)
resolved = validate_contract(root, contract)
receipt = {
    "schema_version": "masld-bench-gse244832-label-free-atac-exchange-validation-v1",
    "dataset_id": "gse244832",
    "donors": 18,
    "biological_unit": "donor",
    "condition_values_read": False,
    "rna_assay_opened": False,
    "cross_assay_join_inferred": False,
    "biological_outcome_arrays_read": False,
    "fragment_geometry": contract["window_exchange"]["query_fragment_cut_sites"],
    "family_eligibility_rows": len(resolved["family"]),
    "query_materialization_authorized": False,
    "biological_scoring_authorized": False,
    "champion_claim_allowed": False,
    "status": "passed",
}
write_json_exclusive(target / "validation_receipt.json", receipt)
freeze_tree(
    target,
    {
        "artifact_class": "gse244832_label_free_atac_exchange_validation",
        "dataset_id": "gse244832",
        "condition_values_read": False,
        "cross_assay_join_inferred": False,
        "biological_outcome_arrays_read": False,
        "champion_claim_allowed": False,
        "status": "passed",
    },
)
print(json.dumps(receipt, sort_keys=True))
PY

printf 'GSE244832 label-free exchange validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
