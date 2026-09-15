#!/usr/bin/env bash
#SBATCH --job-name=model-check-230
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
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
CONVERSION="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-21064945"
CONFIG="${ROOT}/config/gse281364_mpralegnet_provenance_audit.toml"
REUSE="${ROOT}/executions/model-check-226-21091805.failed"
TARGET="${ROOT}/executions/model-check-230-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed terminal audit: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite terminal audit: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"

for job in "21091805 FAILED 1:0" "21092032 FAILED 1:0"; do
  read -r job_id expected_state expected_exit <<< "${job}"
  state="$(sacct -X -j "${job_id}" --format=State -n -P | sed -n '1p' | cut -d'|' -f1)"
  exit_code="$(sacct -X -j "${job_id}" --format=ExitCode -n -P | sed -n '1p' | cut -d'|' -f1)"
  if [[ "${state}" != "${expected_state}" || "${exit_code}" != "${expected_exit}" ]]; then
    printf 'Preserved failure disposition differs: %s\n' "${job_id}" >&2
    exit 1
  fi
done

for binding in \
  "6f77ec23b4990efe416a3f69fcacb9ccce942b2fb062f17855fb9bd32e26ca7f executions/mpralegnet-admission-21066221" \
  "33c3ea9fe134157b9f8f29396cd9d0d565243a9f6e44db3f458cd970c486d32c executions/gse281364-mpralegnet-fixture-21069059" \
  "a59c82c5093a61025f39d5fb12c61820711f58e001ee9bc7858bd505140ec5f9 executions/model-work-102-21069877" \
  "fdc1db620152757f74c07904487a62381fcd182a532669b93b0aa1ca5cdfc2ac executions/model-check-225-21091492"; do
  read -r expected relative <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${ROOT}/${relative}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
printf '%s  %s\n' \
  "b2996006b72c41b4232d0c2b01d445c5ad040760a7e32b5e29ce1fe4d0cefb24" \
  "${ROOT}/executions/model-check-229-21092032.failed/validation/audit.stderr.txt" \
  | sha256sum --check --strict
printf '%s  %s\n' \
  "9acaa26ffeb3fa06c587426768a75215775625be1d5b869bc390f88a8b6955bf" \
  "${ROOT}/executions/model-check-229-21092032.failed/validation/tests.stderr.txt" \
  | sha256sum --check --strict
printf '%s  %s\n' \
  "34e1c4454341575beeb1158fb7ce08e4258d0b6041f80209498b1fe40a579645" \
  "${REUSE}/downloads/zenodo_record.json" | sha256sum --check --strict
printf '%s  %s\n' \
  "05118d0edf404a43aeca5e5ca87ec7b3" \
  "${REUSE}/downloads/final_dump.zip" | md5sum --check --strict
printf '%s  %s\n' \
  "de01b24b9cc5960a763df07f0358a19b" \
  "${REUSE}/downloads/human_legnet-main.zip" | md5sum --check --strict

module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages:${ROOT}/src:${ROOT}"
python -m unittest -v tests.unit.test_audit_mpralegnet_checkpoint_provenance \
  > "${STAGE}/validation/tests.stdout.txt" \
  2> "${STAGE}/validation/tests.stderr.txt"
python "${ROOT}/scripts/audit_mpralegnet_checkpoint_provenance.py" \
  --project-root "${ROOT}" \
  --config "${CONFIG}" \
  --zenodo-record "${REUSE}/downloads/zenodo_record.json" \
  --final-dump "${REUSE}/downloads/final_dump.zip" \
  --code-archive "${REUSE}/downloads/human_legnet-main.zip" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json" \
  2> "${STAGE}/validation/audit.stderr.txt"
python -m pip freeze --all > "${STAGE}/validation/pip-freeze.txt"
python --version > "${STAGE}/validation/python.version.txt" 2>&1
sha256sum \
  "${ROOT}/config/gse281364_mpralegnet_provenance_audit.toml" \
  "${ROOT}/scripts/audit_mpralegnet_checkpoint_provenance.py" \
  "${ROOT}/tests/unit/test_audit_mpralegnet_checkpoint_provenance.py" \
  "${ROOT}/slurm/freeze_mpralegnet_terminal_provenance_cpu.sbatch" \
  > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "audit/provenance/receipt.json").read_text())
inventory = json.loads(
    (root / "audit/provenance/restricted_serialization_inventory.json").read_text()
)
if (
    receipt["status"]
    != "blocked_official_weights_only_scanner_unsupported_POP_opcode"
    or receipt["provenance_equivalence_pass"]
    or receipt["descriptive_native_output_scoring_eligible"]
    or receipt["checkpoint_state_tensors_bitwise_equal"] is not None
    or receipt["torch_load_invoked"]
    or receipt["weights_only_load_attempted"]
    or receipt["weights_only_false_used"]
    or receipt["custom_deserializer_used"]
    or receipt["unrestricted_deserializer_used"]
    or receipt["outcomes_read"]
    or receipt["metrics_calculated"]
    or receipt["model_fit"]
    or receipt["predictions_generated"]
    or inventory["status"] != receipt["status"]
    or (root / "audit/contract").exists()
    or (root / "audit/provenance/checkpoint_tensor_equivalence.tsv").exists()
):
    raise SystemExit("terminal provenance receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "mpralegnet_terminal_checkpoint_provenance_audit",
        "dataset_id": "gse281364",
        "model_id": "mpralegnet_hepg2_test1_val2",
        "terminal_reason": receipt["status"],
        "provenance_equivalence_pass": False,
        "pre_scoring_contract_frozen": False,
        "torch_load_invoked": False,
        "custom_deserializer_used": False,
        "unrestricted_deserializer_used": False,
        "outcomes_read": False,
        "metrics_calculated": False,
        "model_fit": False,
        "predictions_generated": False,
        "standalone_champion_eligible": False,
        "status": "terminal_blocked",
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'MPRALegNet terminal provenance audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
