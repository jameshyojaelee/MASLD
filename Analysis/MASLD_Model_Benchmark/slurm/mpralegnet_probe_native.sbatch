#!/usr/bin/env bash
#SBATCH --job-name=masld-mpralegnet-native
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
export CUBLAS_WORKSPACE_CONFIG=:4096:8

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
BASE_SHA=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_SHA=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
CONVERSION="${ROOT}/executions/environments/dna-lm-torch2.6-conversion-21064945"
CONVERSION_SHA=5fed22c975200495ea53b38b09f5286cf13a9230a138cb72a0c758c85d5cb9d4
ADMISSION="${ROOT}/executions/mpralegnet-admission-21066221"
ADMISSION_SHA=6f77ec23b4990efe416a3f69fcacb9ccce942b2fb062f17855fb9bd32e26ca7f
FIXTURE="${ROOT}/executions/dna-lm-sequence-fixture-21064549"
FIXTURE_SHA=c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2
TARGET="${ROOT}/executions/mpralegnet-native-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed MPRALegNet native probe: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite MPRALegNet native probe: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/offline-guard"
for binding in \
  "${BASE_SHA} ${BASE}" "${TORCH_SHA} ${TORCH}" \
  "${CONVERSION_SHA} ${CONVERSION}" "${ADMISSION_SHA} ${ADMISSION}" \
  "${FIXTURE_SHA} ${FIXTURE}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${BASE}" "${TORCH}" "${CONVERSION}" "${ADMISSION}" "${FIXTURE}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("MPRALegNet native frozen inputs verified")
PY
module purge

cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket
_socket = socket.socket
class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for MPRALegNet execution")
        return super().__new__(cls, family, *args, **kwargs)
socket.socket = GuardedSocket
def blocked_connection(*args, **kwargs):
    raise RuntimeError("network disabled for MPRALegNet execution")
socket.create_connection = blocked_connection
PY
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${STAGE}/offline-guard:${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages"
python "${ROOT}/scripts/mpralegnet_probe_native.py" \
  --model-source "${ADMISSION}/sources/author_model.py" \
  --config "${ADMISSION}/metadata/config.json" \
  --checkpoint "${ADMISSION}/checkpoint/hepg2_best_model_test1_val2.safetensors" \
  --manifest "${FIXTURE}/fixture/sequence_manifest.tsv" \
  --fasta "${FIXTURE}/fixture/common_6000.alleles.fa.gz" \
  --output "${STAGE}/predictions" \
  > "${STAGE}/validation/probe.stdout.txt" \
  2> "${STAGE}/validation/probe.stderr.txt"
python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/mpralegnet_probe_native.py" \
  "${ROOT}/slurm/mpralegnet_probe_native.sbatch" \
  "${ADMISSION}/ARTIFACTS.json" "${FIXTURE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "predictions/probe_receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["checkpoint_sha256"] != "470dc7bfd3f0912c91f307a5f4019598072b8786294c46fc501b826c030cbe39"
    or receipt["parameter_and_buffer_elements"] != 1_330_548
    or receipt["input_length_bp"] != 230
    or receipt["channel_order"] != "A_G_C_T"
    or receipt["forward_variant_index0"] != 115
    or receipt["reverse_complement_variant_index0"] != 114
    or receipt["allele_effect_sign"] != "ALT_minus_REF"
    or not receipt["strict_checkpoint_restore"]
    or not receipt["batchnorm_running_state_unchanged"]
    or receipt["deterministic_repeat_max_abs_diff"] != 0.0
    or not receipt["input_gradient_backward_finite"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_outcomes_loaded"]
    or receipt["standalone_champion_eligible"]
):
    raise SystemExit("MPRALegNet native probe receipt differs")
freeze_tree(root, {
    "artifact_class": "mpralegnet_native_runtime_probe",
    "admission_artifacts_sha256": "6f77ec23b4990efe416a3f69fcacb9ccce942b2fb062f17855fb9bd32e26ca7f",
    "fixture_artifacts_sha256": "c5018dfd39f4ebde636e76ae4d396bfb3c0b326a8f84894c0b4115f7e8b847b2",
    "model_id": "mpralegnet_hepg2_test1_val2",
    "checkpoint_provenance": receipt["checkpoint_provenance"],
    "input_contract": "230bp_AGCT_variant_index115_RC_mean_ALT_minus_REF",
    "observed_outcomes_loaded": False,
    "sealed_outcomes_loaded": False,
    "standalone_champion_eligible": False,
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'MPRALegNet native probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
