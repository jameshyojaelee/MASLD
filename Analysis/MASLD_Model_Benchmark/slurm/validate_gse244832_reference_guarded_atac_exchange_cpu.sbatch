#!/usr/bin/env bash
#SBATCH --job-name=model-check-334
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
TARGET="${ROOT}/executions/model-check-334-${SLURM_JOB_ID}"
if [[ -e "${TARGET}" ]]; then
  printf 'Refusing to overwrite reference-guard validation: %s\n' "${TARGET}" >&2
  exit 1
fi
mkdir -p "${TARGET}"
cd "${ROOT}"

for pair in \
  "402ade8ffd26c450f4a763591870868620692f96f97c92310dac9174570d02b5 config/evaluation/gse244832_label_free_atac_exchange.toml" \
  "da94da73720a1e799e49edcfc261f37ef53cf6bcbb6abd4c92773d3e84fef368 config/evaluation/gse244832_atac_reference_guard_20260825.json" \
  "881391442717d4eccc71a3aa0e1861fa8d1c027f19796a992da660eed17ff3c2 scripts/build_gse244832_reference_guarded_atac_exchange.py" \
  "5ce29f4faa7735c143c1b5f0244bdf5a01a7ccd986f9c26013964c61c6992fe3 tests/unit/test_build_gse244832_reference_guarded_atac_exchange.py"; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" -m unittest -v \
  tests.unit.test_build_gse244832_reference_guarded_atac_exchange \
  tests.unit.test_gse244832_label_free_atac_exchange
"${PYTHON}" -m pip freeze > "${TARGET}/pip-freeze.txt"

PYTHONPATH="${ROOT}/src:${ROOT}" "${PYTHON}" - "${ROOT}" "${TARGET}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, write_json_exclusive
from scripts.build_gse244832_reference_guarded_atac_exchange import validate_guard

root = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2]).resolve(strict=True)
contract = root / "config/evaluation/gse244832_label_free_atac_exchange.toml"
guard = root / "config/evaluation/gse244832_atac_reference_guard_20260825.json"
builder = root / "scripts/build_gse244832_reference_guarded_atac_exchange.py"
test = root / "tests/unit/test_build_gse244832_reference_guarded_atac_exchange.py"
_contract, _guard, resolved = validate_guard(root, contract, guard, builder, test)
receipt = {
    "schema_version": "masld-bench-gse244832-reference-guard-validation-v1",
    "dataset_id": "gse244832",
    "family_rows": len(resolved["family"]),
    "non_observed_atac_executable_rows": sum(
        row["input_regime"] != "observed_atac"
        and not row["eligibility"].startswith("ineligible_")
        for row in resolved["family"]
    ),
    "sequence_execution_authorized": False,
    "observed_atac_query_materialization_authorized": True,
    "condition_values_read": False,
    "rna_assay_opened": False,
    "biological_outcome_arrays_read": False,
    "champion_claim_allowed": False,
    "status": "passed",
}
if receipt["non_observed_atac_executable_rows"] != 0:
    raise SystemExit("reference guard left non-observed-ATAC family executable")
write_json_exclusive(target / "validation_receipt.json", receipt)
freeze_tree(
    target,
    {
        "artifact_class": "gse244832_reference_guarded_atac_exchange_validation",
        "dataset_id": "gse244832",
        "sequence_execution_authorized": False,
        "observed_atac_query_materialization_authorized": True,
        "condition_values_read": False,
        "rna_assay_opened": False,
        "biological_outcome_arrays_read": False,
        "champion_claim_allowed": False,
        "status": "passed",
    },
)
print(json.dumps(receipt, sort_keys=True))
PY

printf 'GSE244832 reference-guard validation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
