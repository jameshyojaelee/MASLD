#!/usr/bin/env bash
#SBATCH --job-name=masld-sei-native
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
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
SOURCE="${ROOT}/executions/alphagenome-sei-admission-sources-21064512"
SOURCE_SHA=84ae3470dca123c20904e315ef9926d5a139fb9fe6d99c51f6aaf694ef6e6785
SCORING="${ROOT}/executions/sei-scoring-sources-21065483"
SCORING_SHA=529a2034a4a8168d8e436c7044364c65735bb274a43b9652d65000110562c5f7
WEIGHTS="${ROOT}/executions/sei-weights-only-conversion-21065468"
WEIGHTS_SHA=2b19a580f496db0a30a6174cc4f7fe65b0cf399e5dde94c8087665d374ba892b
FIXTURE="${ROOT}/executions/alphagenome-sei-sequence-fixture-21064515"
FIXTURE_SHA=e87468da0f64da1c022751834b624471221aec3ba03fcd17d6dbf756e3b8472e
TARGET="${ROOT}/executions/sei-native-probe-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Sei native stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Sei native artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation" "${STAGE}/offline-guard"
for binding in \
  "${BASE_SHA} ${BASE}" \
  "${TORCH_SHA} ${TORCH}" \
  "${CONVERSION_SHA} ${CONVERSION}" \
  "${SOURCE_SHA} ${SOURCE}" \
  "${SCORING_SHA} ${SCORING}" \
  "${WEIGHTS_SHA} ${WEIGHTS}" \
  "${FIXTURE_SHA} ${FIXTURE}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" \
    | sha256sum --check --strict
done

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - \
  "${BASE}" "${TORCH}" "${CONVERSION}" "${SOURCE}" "${SCORING}" "${WEIGHTS}" "${FIXTURE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("Sei native frozen inputs verified")
PY
module purge

cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket

_socket = socket.socket

class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for checkpoint execution")
        return super().__new__(cls, family, *args, **kwargs)

socket.socket = GuardedSocket

def blocked_connection(*args, **kwargs):
    raise RuntimeError("network disabled for checkpoint execution")

socket.create_connection = blocked_connection
PY
module load python/3.12.3-GCCcore-13.3.0
export PYTHONPATH="${STAGE}/offline-guard:${CONVERSION}/site-packages:${TORCH}/site-packages:${BASE}/site-packages"
python "${ROOT}/scripts/sei_probe_native.py" \
  --model-source "${SOURCE}/sources/sei_model_py" \
  --scoring-source "${SCORING}/sources/utils.py" \
  --checkpoint "${WEIGHTS}/converted/sei.model.safetensors" \
  --projection "${WEIGHTS}/native/projvec_targets.npy" \
  --histone-indices "${WEIGHTS}/native/histone_inds.npy" \
  --manifest "${FIXTURE}/fixture/sequence_manifest.tsv" \
  --fasta "${FIXTURE}/fixture/sei.alleles.fa.gz" \
  --output "${STAGE}/predictions" \
  > "${STAGE}/validation/probe.stdout.txt" \
  2> "${STAGE}/validation/probe.stderr.txt"
python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/sei_probe_native.py" \
  "${ROOT}/slurm/sei_probe_native.sbatch" \
  "${SOURCE}/ARTIFACTS.json" \
  "${SCORING}/ARTIFACTS.json" \
  "${WEIGHTS}/ARTIFACTS.json" \
  "${FIXTURE}/ARTIFACTS.json" \
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
    or receipt["parameter_count"] != 889_979_983
    or receipt["profile_shape"] != [12, 21907]
    or receipt["native_variant_score_shape"] != [3, 40]
    or receipt["deterministic_repeat_max_abs_diff"] != 0.0
    or not receipt["input_gradient_backward_finite"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_outcomes_loaded"]
    or receipt["champion_eligible"]
):
    raise SystemExit("Sei native probe receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "sei_native_runtime_probe",
        "source_artifacts_sha256": "84ae3470dca123c20904e315ef9926d5a139fb9fe6d99c51f6aaf694ef6e6785",
        "scoring_artifacts_sha256": "529a2034a4a8168d8e436c7044364c65735bb274a43b9652d65000110562c5f7",
        "weights_artifacts_sha256": "2b19a580f496db0a30a6174cc4f7fe65b0cf399e5dde94c8087665d374ba892b",
        "fixture_artifacts_sha256": "e87468da0f64da1c022751834b624471221aec3ba03fcd17d6dbf756e3b8472e",
        "parameter_count": 889_979_983,
        "profile_output_count": 21_907,
        "native_variant_output_count": 40,
        "terms": "academic_and_research_use_only_restricted_comparator",
        "observed_outcomes_loaded": False,
        "sealed_outcomes_loaded": False,
        "champion_eligible": False,
    },
)
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Sei native runtime probe: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
