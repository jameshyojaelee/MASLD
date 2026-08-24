#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-schema
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ACQUISITION="${ROOT}/executions/upstream-acquisitions/corgi-d84c000-zenodo18630048-21064463"
ACQUISITION_ARTIFACTS_SHA256=73aad9ab9ad3e080a2ccaa0cee1a1d727d126a454f0e69f74f343a900b168c0f
INVENTORY_SHA256=2e92448584c148b7e8893245c30e2ffee6c956ca571c908fb4b7279251ba314f
PYTHON=/gpfs/commons/home/jameslee/scImmuneATLAS/.micromamba/envs/immune-atlas/bin/python
TARGET="${ROOT}/executions/corgi-official-checkpoint-schema-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Corgi checkpoint-schema artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/schemas" "${STAGE}/validation"
printf '%s  %s\n' "${ACQUISITION_ARTIFACTS_SHA256}" \
  "${ACQUISITION}/ARTIFACTS.json" | sha256sum --check --strict
printf '%s  %s\n' "${INVENTORY_SHA256}" \
  "${ROOT}/scripts/corgi_inventory_checkpoint.py" | sha256sum --check --strict

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${ACQUISITION}" \
  "${STAGE}/validation/acquisition-binding.json" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import verify_frozen_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
manifest = verify_frozen_tree(root)
source = root / "source/corgi-d84c00082da3d04eb05169e290e53f1c7da2b9c5"
predict = (source / "corgi/predict.py").read_text(encoding="utf-8")
if (
    "class corgi_pretrained" not in predict
    or "class corgiplus_pretrained" not in predict
    or "raise NotImplementedError" not in predict
):
    raise SystemExit("corrected Corgi public-API source audit differs")
payload = {
    "artifact_manifest_sha256": sha256_file(root / "ARTIFACTS.json"),
    "artifact_metadata": manifest["metadata"],
    "regular_high_level_api_present": True,
    "corgi_plus_placeholder_present": True,
    "released_loader_weights_only": False,
    "released_loader_strict": False,
    "upstream_acquisition_source_contract_caveat": (
        "The frozen acquisition source-contract detector searched for def rather "
        "than class and retained pre-publication staging path strings. This bound "
        "audit corrects those two metadata fields; source and checkpoint bytes are "
        "unchanged and recursively verified."
    ),
    "source_sha256": {
        "corgi/predict.py": sha256_file(source / "corgi/predict.py"),
        "corgi/model.py": sha256_file(source / "corgi/model.py"),
        "corgi/config.py": sha256_file(source / "corgi/config.py"),
        "corgi/config_corgiplus.py": sha256_file(source / "corgi/config_corgiplus.py"),
    },
}
Path(sys.argv[2]).write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY

module purge
"${PYTHON}" "${ROOT}/scripts/corgi_inventory_checkpoint.py" \
  --checkpoint "${ACQUISITION}/downloads/corgi_model.pt" \
  --expected-size 2350451066 \
  --expected-md5 b10ddca8bb32eae3b67116242cfe1871 \
  --expected-sha256 cd51539f01de1e66a3a4f9b89e772bdeae07df2f0178d5b692f2368db2f4d2a4 \
  --output "${STAGE}/schemas/corgi_regular.json"
"${PYTHON}" "${ROOT}/scripts/corgi_inventory_checkpoint.py" \
  --checkpoint "${ACQUISITION}/downloads/corgiplus_model.pt" \
  --expected-size 2384030158 \
  --expected-md5 f7dbe2ca301fe64814bebcd4fc3377bd \
  --expected-sha256 01855614bfeffa72d5bcaaa28729c022ccfbc8f5fdcc0add8e04ba4533e94ccd \
  --output "${STAGE}/schemas/corgi_plus.json"
"${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

root = Path(sys.argv[1])
regular = json.loads((root / "schemas/corgi_regular.json").read_text(encoding="utf-8"))
plus = json.loads((root / "schemas/corgi_plus.json").read_text(encoding="utf-8"))
if regular["wrapper_path"] != ["model_state_dict"]:
    raise SystemExit("Regular Corgi checkpoint wrapper differs")
if plus["wrapper_path"] != ["model_state_dict"]:
    raise SystemExit("Corgi+ checkpoint wrapper differs")
if regular["auxiliary_key_records"]:
    raise SystemExit("Regular Corgi unexpectedly contains auxiliary parameters")
aux_input = [
    record for record in plus["state"]
    if str(record["key"]).endswith("aux_encoder.0.weight")
]
if len(aux_input) != 1 or len(aux_input[0]["shape"]) != 3:
    raise SystemExit("Corgi+ auxiliary input tensor is not uniquely identifiable")
aux_width = int(aux_input[0]["shape"][1])
if aux_width not in {22, 24, 26}:
    raise SystemExit(f"Corgi+ auxiliary width differs: {aux_width}")
summary = {
    "schema_version": "masld-bench-corgi-official-checkpoint-schema-v1",
    "status": "pass",
    "safe_loader": "torch.load(weights_only=True,map_location=cpu,mmap=True)",
    "regular": {
        "sha256": regular["checkpoint_sha256"],
        "state_key_count": regular["state_key_count"],
        "total_numel": regular["total_numel"],
        "total_storage_bytes": regular["total_storage_bytes"],
        "auxiliary_key_count": 0,
    },
    "corgi_plus": {
        "sha256": plus["checkpoint_sha256"],
        "state_key_count": plus["state_key_count"],
        "total_numel": plus["total_numel"],
        "total_storage_bytes": plus["total_storage_bytes"],
        "auxiliary_input_width": aux_width,
        "auxiliary_input_tensor": aux_input[0],
    },
    "tensor_values_exported": False,
    "corgi_source_imported": False,
    "strict_architecture_restore_performed": False,
    "native_numeric_parity_performed": False,
    "evaluator_outcomes_exposed": False,
    "regular_admission_status": "BLOCKED_PENDING_LICENSE_STRICT_RESTORE_RUNTIME_REFERENCE_AND_NUMERIC_PARITY",
    "corgi_plus_admission_status": "BLOCKED_INCOMPLETE_RELEASED_INPUT_AND_NUMERIC_PARITY",
}
(root / "summary.json").write_text(
    json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(summary, sort_keys=True))
PY
"${PYTHON}" -m pip freeze > "${STAGE}/validation/pip-freeze.txt"
"${PYTHON}" - "${STAGE}/validation/runtime.json" <<'PY'
from pathlib import Path
import json
import platform
import sys
import torch

Path(sys.argv[1]).write_text(json.dumps({
    "python": platform.python_version(),
    "torch": torch.__version__,
    "torch_cuda": torch.version.cuda,
    "cuda_available": torch.cuda.is_available(),
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

sha256sum \
  "${ROOT}/scripts/corgi_inventory_checkpoint.py" \
  "${ROOT}/slurm/corgi_probe_official_checkpoint_schema.sbatch" \
  > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${ACQUISITION_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

freeze_tree(
    Path(sys.argv[1]),
    metadata={
        "artifact_class": "corgi_official_checkpoint_schema",
        "status": "pass",
        "acquisition_artifacts_sha256": sys.argv[2],
        "safe_weights_only_load": True,
        "strict_architecture_restore_performed": False,
        "native_numeric_parity_performed": False,
        "evaluator_outcomes_exposed": False,
    },
)
PY
mv "${STAGE}" "${TARGET}"
printf 'Corgi official checkpoint schema: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
