#!/usr/bin/env bash
#SBATCH --job-name=masld-epibert-exec-runtime
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE="${ROOT}/executions/environments/epibert-tensorflow-2.13.0-21064371"
BASE_ARTIFACTS_SHA256=10a8691e9468dfa44a25e55770b477e77d2ebed12839d49203f1a386596af27f
IMAGE="${BASE}/tensorflow-tensorflow-2.13.0-gpu.sif"
SOURCES="${ROOT}/executions/epibert-epcotv2-admission-sources-21064417"
SOURCES_ARTIFACTS_SHA256=26e12915b05c2727c45eca091f85ec89f96c205e817c6ec9b64f181a563e7577
TARGET="${ROOT}/executions/environments/epibert-executable-tf2.13-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite EpiBERT executable runtime: %s\n' "${TARGET}" >&2
  exit 1
fi
printf '%s  %s\n' "${BASE_ARTIFACTS_SHA256}" "${BASE}/ARTIFACTS.json" \
  | sha256sum --check --strict
printf '%s  %s\n' "${SOURCES_ARTIFACTS_SHA256}" "${SOURCES}/ARTIFACTS.json" \
  | sha256sum --check --strict
install -d -m 0750 \
  "${STAGE}/source" \
  "${STAGE}/checkpoints/pretrained_model1" \
  "${STAGE}/checkpoints/pretrained_model2" \
  "${STAGE}/checkpoints/fine_tuned_rampage" \
  "${STAGE}/upstream_checkpoint_state" \
  "${STAGE}/wheels" \
  "${STAGE}/site-packages" \
  "${STAGE}/validation"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${BASE}" "${SOURCES}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import verify_frozen_tree

for value in sys.argv[1:]:
    verify_frozen_tree(Path(value))
PY
python - "${SOURCES}" "${STAGE}" <<'PY'
from pathlib import Path, PurePosixPath
import hashlib
import json
import shutil
import sys
import tarfile

source_root = Path(sys.argv[1]).resolve(strict=True)
target = Path(sys.argv[2]).resolve(strict=True)
code_archive = source_root / "downloads/epibert_code.tar.gz"
model_archive = source_root / "downloads/epibert_models.tar.gz"
code_inventory = json.loads(
    (source_root / "inventories/epibert_code.inventory.json").read_text()
)
model_inventory = json.loads(
    (source_root / "inventories/epibert_models.inventory.json").read_text()
)

def safe_members(handle):
    members = handle.getmembers()
    seen = set()
    for member in members:
        path = PurePosixPath(member.name)
        if (
            path.is_absolute()
            or ".." in path.parts
            or member.name in seen
            or not (member.isfile() or member.isdir())
        ):
            raise SystemExit(f"unsafe archive member: {member.name!r}")
        seen.add(member.name)
    return members

with tarfile.open(code_archive, "r:gz") as handle:
    members = safe_members(handle)
    handle.extractall(target / "source", members=members, filter="data")

inventory_by_path = {
    row["path"]: row
    for row in model_inventory["members"]
    if row["type"] == "regular_file"
}

def digest_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

selected = {
    "epibert_models/pretrained/model1/ckpt-45.data-00000-of-00001": target / "checkpoints/pretrained_model1/ckpt-45.data-00000-of-00001",
    "epibert_models/pretrained/model1/ckpt-45.index": target / "checkpoints/pretrained_model1/ckpt-45.index",
    "epibert_models/pretrained/model2/ckpt-45.data-00000-of-00001": target / "checkpoints/pretrained_model2/ckpt-45.data-00000-of-00001",
    "epibert_models/pretrained/model2/ckpt-45.index": target / "checkpoints/pretrained_model2/ckpt-45.index",
    "epibert_models/fine_tuned/ckpt-16.data-00000-of-00001": target / "checkpoints/fine_tuned_rampage/ckpt-16.data-00000-of-00001",
    "epibert_models/fine_tuned/ckpt-16.index": target / "checkpoints/fine_tuned_rampage/ckpt-16.index",
    "epibert_models/pretrained/model1/checkpoint": target / "upstream_checkpoint_state/pretrained_model1.checkpoint",
    "epibert_models/pretrained/model2/checkpoint": target / "upstream_checkpoint_state/pretrained_model2.checkpoint",
    "epibert_models/fine_tuned/checkpoint": target / "upstream_checkpoint_state/fine_tuned_rampage.checkpoint",
}
with tarfile.open(model_archive, "r:gz") as handle:
    members = {member.name: member for member in safe_members(handle)}
    for source, destination in selected.items():
        member = members.get(source)
        expected = inventory_by_path.get(source)
        if member is None or expected is None or member.size != expected["size_bytes"]:
            raise SystemExit(f"selected checkpoint member differs: {source}")
        fileobj = handle.extractfile(member)
        if fileobj is None:
            raise SystemExit(f"selected checkpoint member is unreadable: {source}")
        with destination.open("xb") as writer:
            shutil.copyfileobj(fileobj, writer, length=8 * 1024 * 1024)
        digest = digest_file(destination)
        if digest != expected["sha256"]:
            raise SystemExit(f"selected checkpoint hash differs: {source}")

receipt = {
    "schema_version": "masld-bench-epibert-selected-checkpoint-extraction-v1",
    "status": "pass",
    "source_revision": "dae61b434f885991c29e2fde3fbf7f70d0ff9b80",
    "selected_checkpoint_members": {
        source: {
            "destination": str(destination.relative_to(target)),
            "sha256": inventory_by_path[source]["sha256"],
            "size_bytes": inventory_by_path[source]["size_bytes"],
        }
        for source, destination in selected.items()
    },
    "nested_duplicate_model2_extracted": False,
    "checkpoint_state_files_trusted_for_restore_prefix": False,
    "explicit_restore_prefixes": {
        "pretrained_model1": "checkpoints/pretrained_model1/ckpt-45",
        "pretrained_model2": "checkpoints/pretrained_model2/ckpt-45",
        "fine_tuned_rampage": "checkpoints/fine_tuned_rampage/ckpt-16",
    },
    "reason_checkpoint_state_ignored": "released_checkpoint_state_lists_missing_later_steps; only the explicit released index/data prefixes are admissible",
}
(target / "validation/selected-checkpoint-extraction.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
)
PY

module purge
module load singularity/4.2.2
singularity exec --cleanenv \
  --bind "${STAGE}:${STAGE}:rw" \
  "${IMAGE}" python -m pip download \
    --disable-pip-version-check \
    --only-binary=:all: \
    --no-deps \
    --dest "${STAGE}/wheels" \
    tensorflow-addons==0.21.0 \
    einops==0.6.1 \
    typeguard==2.13.0
sha256sum "${STAGE}"/wheels/* > "${STAGE}/validation/wheels.sha256"
singularity exec --cleanenv \
  --bind "${STAGE}:${STAGE}:rw" \
  "${IMAGE}" python -m pip install \
    --disable-pip-version-check \
    --no-index \
    --no-deps \
    --find-links "${STAGE}/wheels" \
    --target "${STAGE}/site-packages" \
    tensorflow-addons==0.21.0 \
    einops==0.6.1 \
    typeguard==2.13.0

SOURCE="${STAGE}/source/EpiBERT-dae61b434f885991c29e2fde3fbf7f70d0ff9b80"
singularity exec --cleanenv \
  --bind "${STAGE}:${STAGE}:ro" \
  --env "PYTHONPATH=${SOURCE}:${STAGE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python - "${STAGE}" <<'PY' \
  > "${STAGE}/validation/executable-runtime-probe.json"
from pathlib import Path
import importlib
import json
import sys
import tensorflow as tf
import tensorflow_addons as tfa
import einops

root = Path(sys.argv[1]).resolve(strict=True)
atac = importlib.import_module("src.models.epibert_atac_pretrain")
rampage = importlib.import_module("src.models.epibert_rampage_finetune")
checkpoints = {}
for model, relative in {
    "pretrained_model1": "checkpoints/pretrained_model1/ckpt-45",
    "pretrained_model2": "checkpoints/pretrained_model2/ckpt-45",
    "fine_tuned_rampage": "checkpoints/fine_tuned_rampage/ckpt-16",
}.items():
    variables = tf.train.list_variables(str(root / relative))
    if not variables:
        raise SystemExit(f"{model} checkpoint variable list is empty")
    checkpoints[model] = {
        "variable_count": len(variables),
        "total_elements": sum(
            __import__("math").prod(shape) for _name, shape in variables
        ),
        "first_variable": variables[0][0],
        "last_variable": variables[-1][0],
    }
print(json.dumps({
    "status": "pass",
    "python": sys.version.split()[0],
    "tensorflow": tf.__version__,
    "tensorflow_addons": tfa.__version__,
    "einops": einops.__version__,
    "inactive_sonnet_module_imported": False,
    "atac_class": atac.epibert.__name__,
    "rampage_class": rampage.epibert.__name__,
    "checkpoint_inventory": checkpoints,
    "checkpoint_restored": False,
    "model_forward_executed": False,
}, sort_keys=True))
PY
singularity exec --cleanenv \
  --bind "${STAGE}:${STAGE}:ro" \
  --env "PYTHONPATH=${STAGE}/site-packages" \
  --env PYTHONNOUSERSITE=1 \
  "${IMAGE}" python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"

sha256sum \
  "${ROOT}/slurm/epibert_prepare_executable_runtime.sbatch" \
  > "${STAGE}/source.sha256"
module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${BASE_ARTIFACTS_SHA256}" "${SOURCES_ARTIFACTS_SHA256}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation/executable-runtime-probe.json").read_text())
extraction = json.loads(
    (root / "validation/selected-checkpoint-extraction.json").read_text()
)
if (
    probe.get("status") != "pass"
    or probe.get("tensorflow") != "2.13.0"
    or probe.get("tensorflow_addons") != "0.21.0"
    or probe.get("inactive_sonnet_module_imported")
    or probe.get("checkpoint_restored")
    or probe.get("model_forward_executed")
    or extraction.get("nested_duplicate_model2_extracted")
    or extraction.get("checkpoint_state_files_trusted_for_restore_prefix")
):
    raise SystemExit("EpiBERT executable runtime probe differs")
freeze_tree(
    root,
    {
        "artifact_class": "epibert_executable_runtime",
        "base_runtime_artifacts_sha256": sys.argv[2],
        "admission_sources_artifacts_sha256": sys.argv[3],
        "tensorflow": "2.13.0",
        "tensorflow_addons": "0.21.0",
        "code_revision": "dae61b434f885991c29e2fde3fbf7f70d0ff9b80",
        "checkpoint_members": sorted(probe["checkpoint_inventory"]),
        "checkpoint_state_files_trusted_for_restore_prefix": False,
        "checkpoint_restored": False,
        "model_forward_executed": False,
        "project_data_read": False,
        "outcomes_read": False,
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'EpiBERT executable runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
