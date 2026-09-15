#!/usr/bin/env bash
#SBATCH --job-name=model-data-311
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
TARGET=${ROOT}/executions/protein-transport-fixture-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve EXIT
for p in "${TARGET}" "${STAGE}" "${FAILED}"; do
  if [[ -e "${p}" ]]; then printf 'Refusing to overwrite: %s\n' "${p}" >&2; exit 1; fi
done
install -d -m 0750 "${STAGE}" "${STAGE}/audit" "${STAGE}/validation"

cd "${ROOT}"
module purge
module load python/3.12.3-GCCcore-13.3.0

ACTIVATION=${ROOT}/executions/pxd051911-activation-readiness-21109461

# 1. Pin every frozen input by hash before reading a byte of it.
printf '%s  %s\n' \
  "8488568d468f4fe15dff152cd1529f069b69b9ef6f27708a518fdfc757b681a7" "${ACTIVATION}/ARTIFACTS.json" \
  | sha256sum --check --strict

for pair in \
  "b2a5fea26985d7b10fe553afadbb693ea0bc47b6093d6567971cf4d2e5efade7 audit/readiness/participant_join.tsv" \
  "c4fafdfb2339d54d17a91779098012e08448aa5ebcc8a865988dcc49f92c969e source/liver_protein_quant.txt" \
  "074bea52448aa4993fe72e606dde3a44ff4107657f915fbf84728311ee7c0f6d source/scwat_protein_quant.txt" \
  "c957481cccb01685505c39d3c6ffffcc1eeebf3c3f5e1cb5114f6e019861079b source/owat_protein_quant.txt" \
  "340874d4f720f909fcfca428fed6da16ba3f7cfdf4d9fd53dd6e0378198ba2bc source/meta_data.txt" ; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ACTIVATION}/${relative}" | sha256sum --check --strict
done

# The registry files this TaskSpec is additive to must be byte-identical.
for pair in \
  "acd9f333bd4d639d6df6b3017c1bf9a040ee921195b59b01fee934f1726a0495 config/evaluation/family_native_tournament.toml" \
  "64d6cb94505e52a83c6e928df811f4b5261574fe3a81ca60a551fac9807138eb config/datasets/pxd051911.toml" ; do
  read -r expected relative <<< "${pair}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}" | sha256sum --check --strict
done

# 2. Tests. The suite carries known pre-existing failures, so guard set -e and let
#    the baseline check adjudicate whether this campaign introduced a NEW one.
set +e
PYTHONPATH=${ROOT}:${ROOT}/src python -m unittest -v \
  tests.unit.test_build_pxd051911_protein_transport_fixture \
  tests.contract.test_protein_transport_task_contract \
  tests.unit.test_check_contract_suite_baseline \
  tests.contract.test_registry_files \
  > "${STAGE}/validation/unittest.txt" 2>&1
set -e

# 3. No-NEW-requirements-failure check. Baseline is 4 pre-existing failures.
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/check_contract_suite_baseline.py \
  --unittest-output "${STAGE}/validation/unittest.txt" \
  --require-module test_build_pxd051911_protein_transport_fixture \
  --require-module test_protein_transport_task_contract \
  --require-module test_registry_files \
  | tee "${STAGE}/validation/contract_baseline.txt"

# 4. Build the fixture.
PYTHONPATH=${ROOT}:${ROOT}/src python scripts/build_pxd051911_protein_transport_fixture.py \
  --source "${ACTIVATION}/source" \
  --join "${ACTIVATION}/audit/readiness/participant_join.tsv" \
  --output "${STAGE}/audit/fixture" \
  > "${STAGE}/validation/build.stdout.txt" 2> "${STAGE}/validation/build.stderr.txt"

# 5. Carry the frozen design documents alongside the fixture they govern.
cp "${ROOT}/config/evaluation/protein_transport_task.toml"                     "${STAGE}/audit/"
cp "${ROOT}/config/evaluation/protein_transport_promotion_gate.json"           "${STAGE}/audit/"
cp "${ROOT}/config/campaigns/protein_transport_power_position_20260825.json"   "${STAGE}/audit/"

# 6. Environment and fixed source record.
pip freeze --all > "${STAGE}/validation/pip-freeze.txt" 2>&1 || true
sha256sum \
  "${ROOT}/scripts/build_pxd051911_protein_transport_fixture.py" \
  "${ROOT}/tests/unit/test_build_pxd051911_protein_transport_fixture.py" \
  "${ROOT}/tests/contract/test_protein_transport_task_contract.py" \
  "${ROOT}/scripts/check_contract_suite_baseline.py" \
  "${ROOT}/config/evaluation/protein_transport_task.toml" \
  "${ROOT}/config/evaluation/protein_transport_promotion_gate.json" \
  "${ROOT}/config/campaigns/protein_transport_power_position_20260825.json" \
  "${ROOT}/config/campaigns/protein_transport_task_registration_20260825.json" \
  "${ROOT}/config/campaigns/ranking_metric_floor_requirement_20260825.json" \
  "${ROOT}/slurm/build_protein_transport_fixture.sbatch" \
  > "${STAGE}/source.sha256"

# 7. Freeze.
PYTHONPATH=${ROOT}/src python3 - "${STAGE}" <<'PY'
import sys
from pathlib import Path
from masld_bench.artifacts import freeze_tree, verify_frozen_tree
stage = Path(sys.argv[1]).resolve(strict=True)
digest = freeze_tree(stage, metadata={
    "artifact_class": "protein_transport_fixture",
    "dataset_id": "pxd051911",
    "task_id": "protein_transport",
    "unit_of_replication": "participant",
    "status": "fixture_built_development_only",
    "champion_claim_allowed": False,
    "rna_protein_pair_inferred": False,
    "plasma_arm_used": False,
    "models_fit": False,
    "models_scored": False,
    "sealed_outcomes_read": False,
})
verify_frozen_tree(stage)
print("ARTIFACTS.json sha256:", digest)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'FROZEN: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
