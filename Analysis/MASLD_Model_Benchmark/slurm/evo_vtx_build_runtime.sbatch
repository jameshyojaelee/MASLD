#!/usr/bin/env bash
#SBATCH --job-name=masld-evo-vtx-runtime
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/scvi-tools-1.5.0.post1-portable-21064579"
BASE_SHA=30814d4944e933d999e8c02c5b521432786b1985cc8794f9436f2a44af5ad1af
TORCH="${ROOT}/executions/environments/pytorch-2.6.0-cu124-portable-21064747"
TORCH_SHA=d4112dbc9574e619edbaec0f2c1503c57e0271fe17517b2e121d711ea21f71a3
ADMISSION="${ROOT}/executions/evo-vtx-runtime-admission-21066350"
ADMISSION_SHA=ed32535056e742537d46d764d5581eeb006d1b84d8f7d4c5b47b99c2a4e031ba
TARGET="${ROOT}/executions/environments/evo-vtx-1.1.0-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Evo/Vortex runtime build: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Evo/Vortex runtime: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/downloads" "${STAGE}/site-packages" \
  "${STAGE}/code" "${STAGE}/offline-guard" "${STAGE}/validation"
for binding in "${BASE_SHA} ${BASE}" "${TORCH_SHA} ${TORCH}" \
  "${ADMISSION_SHA} ${ADMISSION}"; do
  read -r expected directory <<< "${binding}"
  printf '%s  %s\n' "${expected}" "${directory}/ARTIFACTS.json" \
    | sha256sum --check --strict
done
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BASE}" "${TORCH}" "${ADMISSION}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
print("Evo/Vortex runtime frozen inputs verified")
PY
module purge

curl_flags=(--fail --location --retry 6 --retry-all-errors --connect-timeout 30 \
  --proto '=https' --tlsv1.2)
curl "${curl_flags[@]}" --output "${STAGE}/validation/einops-0.8.1.pypi.json" \
  https://pypi.org/pypi/einops/0.8.1/json
module load python/3.11.5-GCCcore-13.2.0
EINOPS_URL=$(python - "${STAGE}/validation/einops-0.8.1.pypi.json" <<'PY'
import json
import sys
data = json.load(open(sys.argv[1]))
matches = [
    row for row in data["urls"]
    if row["filename"] == "einops-0.8.1-py3-none-any.whl"
]
if (
    len(matches) != 1
    or matches[0]["digests"]["sha256"]
    != "919387eb55330f5757c6bea9165c5ff5cfe63a642682ea788a6d472576d81737"
):
    raise SystemExit("einops 0.8.1 wheel identity differs")
print(matches[0]["url"])
PY
)
curl "${curl_flags[@]}" --output "${STAGE}/downloads/einops-0.8.1-py3-none-any.whl" \
  "${EINOPS_URL}"
printf '%s  %s\n%s  %s\n' \
  919387eb55330f5757c6bea9165c5ff5cfe63a642682ea788a6d472576d81737 "${STAGE}/downloads/einops-0.8.1-py3-none-any.whl" \
  0ff9f1db2f9e81e288150b60fd4fe4832b8b992ac2c6c947271b2036ffeb8299 "${ADMISSION}/downloads/vtx-1.1.0-py3-none-any.whl" \
  | sha256sum --check --strict
module purge

module load python/3.12.3-GCCcore-13.3.0
PYTHONPATH="${TORCH}/site-packages:${BASE}/site-packages" python -m pip install \
  --no-deps --no-index --no-compile --target "${STAGE}/site-packages" \
  "${STAGE}/downloads/einops-0.8.1-py3-none-any.whl" \
  "${ADMISSION}/downloads/vtx-1.1.0-py3-none-any.whl" \
  > "${STAGE}/validation/pip-install.stdout.txt" \
  2> "${STAGE}/validation/pip-install.stderr.txt"
tar -xzf "${ADMISSION}/downloads/evo2-source.tar.gz" \
  --no-same-owner --no-same-permissions --strip-components=1 -C "${STAGE}/code"
cat > "${STAGE}/offline-guard/sitecustomize.py" <<'PY'
import socket
_socket = socket.socket
class GuardedSocket(_socket):
    def __new__(cls, family=-1, *args, **kwargs):
        if family in {socket.AF_INET, socket.AF_INET6}:
            raise RuntimeError("network disabled for Evo/Vortex runtime probe")
        return super().__new__(cls, family, *args, **kwargs)
socket.socket = GuardedSocket
def blocked_connection(*args, **kwargs):
    raise RuntimeError("network disabled for Evo/Vortex runtime probe")
socket.create_connection = blocked_connection
PY
export PYTHONPATH="${STAGE}/offline-guard:${STAGE}/site-packages:${TORCH}/site-packages:${BASE}/site-packages"
python - <<'PY' > "${STAGE}/validation/runtime_probe.json"
import importlib.metadata
import importlib.util
import json
import socket
import torch
from vortex.model.model import StripedHyena
from vortex.model.tokenizer import CharLevelTokenizer
from vortex.model.utils import dotdict

blocked = False
try:
    socket.socket(socket.AF_INET, socket.SOCK_STREAM)
except RuntimeError:
    blocked = True
receipt = {
    "schema_version": "masld-bench-evo-vtx-runtime-probe-v1",
    "status": "pass",
    "vtx": importlib.metadata.version("vtx"),
    "einops": importlib.metadata.version("einops"),
    "torch": torch.__version__,
    "torch_cuda": torch.version.cuda,
    "striped_hyena_class_imported": StripedHyena.__name__,
    "tokenizer_class_imported": CharLevelTokenizer.__name__,
    "dotdict_imported": dotdict.__name__,
    "transformer_engine_available": importlib.util.find_spec("transformer_engine") is not None,
    "python_INET_socket_blocked": blocked,
    "model_instantiated": False,
    "checkpoint_loaded": False,
    "model_forward_executed": False,
    "observed_outcomes_loaded": False,
    "sealed_outcomes_loaded": False,
}
if (
    receipt["vtx"] != "1.1.0"
    or receipt["einops"] != "0.8.1"
    or receipt["torch"] != "2.6.0+cu124"
    or receipt["torch_cuda"] != "12.4"
    or receipt["transformer_engine_available"]
    or not blocked
):
    raise SystemExit(f"Evo/Vortex runtime differs: {receipt!r}")
print(json.dumps(receipt, sort_keys=True, indent=2))
PY
python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/slurm/evo_vtx_build_runtime.sbatch" \
  "${ADMISSION}/ARTIFACTS.json" "${BASE}/ARTIFACTS.json" "${TORCH}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"
module purge

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/runtime_probe.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["vtx"] != "1.1.0"
    or receipt["einops"] != "0.8.1"
    or receipt["transformer_engine_available"]
    or not receipt["python_INET_socket_blocked"]
    or receipt["model_instantiated"]
    or receipt["checkpoint_loaded"]
    or receipt["model_forward_executed"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_outcomes_loaded"]
):
    raise SystemExit("Evo/Vortex runtime receipt differs")
freeze_tree(root, {
    "artifact_class": "evo_vtx_runtime",
    "runtime_admission_artifacts_sha256": "ed32535056e742537d46d764d5581eeb006d1b84d8f7d4c5b47b99c2a4e031ba",
    "vtx": "1.1.0",
    "einops": "0.8.1",
    "torch": receipt["torch"],
    "transformer_engine": "not_installed_7b_pure_bfloat16_path",
    "checkpoint_loaded": False,
    "observed_outcomes_loaded": False,
    "sealed_outcomes_loaded": False,
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Evo/Vortex runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
