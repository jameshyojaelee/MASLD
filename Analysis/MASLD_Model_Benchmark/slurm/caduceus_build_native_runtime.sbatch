#!/usr/bin/env bash
#SBATCH --job-name=masld-caduceus-runtime
#SBATCH --partition=cpu
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
SOURCE="${ROOT}/executions/dna-lm-admission-sources-21064540"
SOURCE_SHA=4e1703b03ff28924904053e842b89079d510a80d2a9ad082de5369adaeea75c8
TARGET="${ROOT}/executions/environments/caduceus-native-torch2.2.0-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
export MAMBA_ROOT_PREFIX
MAMBA_ROOT_PREFIX=$(mktemp -d "${TMPDIR:-/tmp}/masld-caduceus-${SLURM_JOB_ID}.XXXXXX")

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Caduceus runtime stage: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Caduceus runtime artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 \
  "${STAGE}/downloads" \
  "${STAGE}/release-metadata" \
  "${STAGE}/validation"
printf '%s  %s\n' "${SOURCE_SHA}" "${SOURCE}/ARTIFACTS.json" \
  | sha256sum --check --strict

micromamba create --yes --prefix "${STAGE}/env" \
  --channel pytorch \
  --channel nvidia \
  --channel defaults \
  python=3.8.18 \
  pip=23.3.1 \
  mkl=2023.1.0=h213fc3f_46344 \
  cuda-version=12.1 \
  pytorch::pytorch=2.2.0=py3.8_cuda12.1_cudnn8.9.2_0 \
  pytorch-cuda=12.1

MAMBA_API='https://api.github.com/repos/state-spaces/mamba/releases/tags/v1.2.0.post1'
CAUSAL_API='https://api.github.com/repos/Dao-AILab/causal-conv1d/releases/tags/v1.2.0.post2'
curl --fail --location --retry 6 --retry-all-errors \
  --output "${STAGE}/release-metadata/mamba-v1.2.0.post1.json" "${MAMBA_API}"
curl --fail --location --retry 6 --retry-all-errors \
  --output "${STAGE}/release-metadata/causal-conv1d-v1.2.0.post2.json" "${CAUSAL_API}"

MAMBA_WHEEL='mamba_ssm-1.2.0.post1+cu122torch2.2cxx11abiFALSE-cp38-cp38-linux_x86_64.whl'
CAUSAL_WHEEL='causal_conv1d-1.2.0.post2+cu122torch2.2cxx11abiFALSE-cp38-cp38-linux_x86_64.whl'
python3 - \
  "${STAGE}/release-metadata/mamba-v1.2.0.post1.json" "${MAMBA_WHEEL}" \
  "${STAGE}/release-metadata/causal-conv1d-v1.2.0.post2.json" "${CAUSAL_WHEEL}" \
  "${STAGE}/validation/release_assets.json" <<'PY'
import json
from pathlib import Path
import sys

pairs = [(Path(sys.argv[1]), sys.argv[2]), (Path(sys.argv[3]), sys.argv[4])]
records = []
for metadata_path, expected_name in pairs:
    metadata = json.loads(metadata_path.read_text())
    matches = [asset for asset in metadata.get("assets", []) if asset.get("name") == expected_name]
    if len(matches) != 1 or metadata.get("draft") or metadata.get("prerelease"):
        raise SystemExit("Caduceus dependency release asset differs")
    asset = matches[0]
    records.append(
        {
            "release_id": metadata["id"],
            "tag_name": metadata["tag_name"],
            "asset_id": asset["id"],
            "name": asset["name"],
            "size": asset["size"],
            "url": asset["browser_download_url"],
            "api_digest": asset.get("digest"),
        }
    )
Path(sys.argv[5]).write_text(json.dumps(records, sort_keys=True, indent=2) + "\n")
print("Caduceus dependency release assets validated")
PY

python3 - "${STAGE}/validation/release_assets.json" "${STAGE}/downloads" <<'PY'
import json
from pathlib import Path
import subprocess
import sys

records = json.loads(Path(sys.argv[1]).read_text())
destination = Path(sys.argv[2])
for record in records:
    subprocess.run(
        [
            "curl",
            "--fail",
            "--location",
            "--retry",
            "6",
            "--retry-all-errors",
            "--output",
            str(destination / record["name"]),
            record["url"],
        ],
        check=True,
    )
    if (destination / record["name"]).stat().st_size != record["size"]:
        raise SystemExit("Caduceus dependency wheel size differs")
PY

cat > "${STAGE}/requirements.lock" <<'EOF'
certifi==2024.2.2
charset-normalizer==3.3.2
einops==0.7.0
huggingface-hub==0.24.7
idna==3.6
numpy==1.24.4
packaging==23.2
PyYAML==6.0.1
regex==2023.12.25
requests==2.31.0
safetensors==0.4.2
tokenizers==0.15.2
tqdm==4.66.1
transformers==4.38.1
urllib3==2.2.1
EOF
"${STAGE}/env/bin/python" -m pip download \
  --disable-pip-version-check \
  --only-binary=:all: \
  --no-deps \
  --dest "${STAGE}/downloads" \
  --requirement "${STAGE}/requirements.lock"
"${STAGE}/env/bin/python" -m pip install \
  --disable-pip-version-check \
  --no-index \
  --no-deps \
  "${STAGE}"/downloads/*.whl

sha256sum "${STAGE}"/downloads/*.whl > "${STAGE}/validation/wheels.sha256"
export HF_HOME="${STAGE}/hf-home"
export HF_HUB_OFFLINE=1
export TRANSFORMERS_OFFLINE=1
install -d -m 0750 "${HF_HOME}"
"${STAGE}/env/bin/python" - <<'PY' > "${STAGE}/validation/import_probe.json"
import json
from importlib.metadata import version
import sys

import einops
import safetensors
import tokenizers
import torch
import transformers

expected = {
    "einops": "0.7.0",
    "huggingface-hub": "0.24.7",
    "mamba-ssm": "1.2.0.post1",
    "causal-conv1d": "1.2.0.post2",
    "numpy": "1.24.4",
    "regex": "2023.12.25",
    "safetensors": "0.4.2",
    "tokenizers": "0.15.2",
    "torch": "2.2.0",
    "transformers": "4.38.1",
}
observed = {name: version(name) for name in expected}
if (
    observed != expected
    or torch.version.cuda != "12.1"
    or torch._C._GLIBCXX_USE_CXX11_ABI
):
    raise SystemExit(f"Caduceus runtime versions differ: {observed!r}")
print(
    json.dumps(
        {
            "python": sys.version,
            "versions": observed,
            "torch_cuda": torch.version.cuda,
            "torch_cxx11_abi": torch._C._GLIBCXX_USE_CXX11_ABI,
            "compiled_mamba_imported": False,
            "checkpoint_loaded": False,
            "model_forward_executed": False,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
)
PY
"${STAGE}/env/bin/python" - "${STAGE}/validation/dependency_metadata.json" <<'PY'
from hashlib import sha256
from importlib.metadata import distribution
import json
from pathlib import Path
import sys

records = []
for name in (
    "certifi",
    "charset-normalizer",
    "einops",
    "huggingface-hub",
    "idna",
    "packaging",
    "PyYAML",
    "mamba-ssm",
    "causal-conv1d",
    "numpy",
    "regex",
    "requests",
    "safetensors",
    "tokenizers",
    "torch",
    "tqdm",
    "transformers",
    "urllib3",
):
    package = distribution(name)
    license_files = []
    for entry in package.files or ():
        value = str(entry)
        if "license" not in value.lower() and "copying" not in value.lower():
            continue
        path = Path(package.locate_file(entry))
        if path.is_file() and not path.is_symlink():
            license_files.append(
                {
                    "path": value,
                    "sha256": sha256(path.read_bytes()).hexdigest(),
                    "size_bytes": path.stat().st_size,
                }
            )
    records.append(
        {
            "name": name,
            "version": package.version,
            "license_metadata": package.metadata.get("License"),
            "license_classifiers": sorted(
                value
                for value in package.metadata.get_all("Classifier", [])
                if value.startswith("License ::")
            ),
            "license_files": license_files,
        }
    )
Path(sys.argv[1]).write_text(json.dumps(records, sort_keys=True, indent=2) + "\n")
PY
"${STAGE}/env/bin/python" - "${STAGE}/env/conda-meta" \
  "${STAGE}/validation/conda_packages.json" <<'PY'
from hashlib import sha256
import json
from pathlib import Path
import sys

source = Path(sys.argv[1])
records = []
for path in sorted(source.glob("*.json")):
    metadata = json.loads(path.read_text())
    records.append(
        {
            "name": metadata.get("name"),
            "version": metadata.get("version"),
            "build": metadata.get("build"),
            "channel": metadata.get("channel"),
            "license": metadata.get("license"),
            "metadata_sha256": sha256(path.read_bytes()).hexdigest(),
        }
    )
if not records:
    raise SystemExit("Caduceus conda package census is empty")
Path(sys.argv[2]).write_text(json.dumps(records, sort_keys=True, indent=2) + "\n")
PY
micromamba list --prefix "${STAGE}/env" --explicit \
  > "${STAGE}/validation/conda-explicit.txt"
"${STAGE}/env/bin/python" -m pip freeze \
  > "${STAGE}/validation/effective-pip-freeze.txt"
micromamba --version > "${STAGE}/validation/micromamba.version.txt"
sha256sum \
  "${ROOT}/slurm/caduceus_build_native_runtime.sbatch" \
  "${SOURCE}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

# Conda environments contain required symlinks, while read-only output-file trees
# forbid live symlinks. Store the complete environment as one checksummed
# payload and prove it still imports after extraction at a different prefix.
install -d -m 0750 "${STAGE}/payload"
tar \
  --sort=name \
  --mtime='UTC 2026-08-24' \
  --owner=0 \
  --group=0 \
  --numeric-owner \
  -C "${STAGE}/env" \
  -czf "${STAGE}/payload/caduceus-native-env.tar.gz" \
  .
RELOCATED=$(mktemp -d "${TMPDIR:-/tmp}/masld-caduceus-relocated-${SLURM_JOB_ID}.XXXXXX")
tar -xzf "${STAGE}/payload/caduceus-native-env.tar.gz" -C "${RELOCATED}"
PYTHONNOUSERSITE=1 "${RELOCATED}/bin/python" - <<'PY' \
  > "${STAGE}/validation/relocation_probe.json"
import json
from importlib.metadata import version
import sys

import torch
import transformers

if (
    version("torch") != "2.2.0"
    or version("transformers") != "4.38.1"
    or torch.version.cuda != "12.1"
    or torch._C._GLIBCXX_USE_CXX11_ABI
):
    raise SystemExit("relocated Caduceus runtime versions differ")
print(
    json.dumps(
        {
            "python": sys.version,
            "prefix_relocated": True,
            "torch": version("torch"),
            "torch_cuda": torch.version.cuda,
            "torch_cxx11_abi": torch._C._GLIBCXX_USE_CXX11_ABI,
            "transformers": version("transformers"),
            "compiled_mamba_imported": False,
            "checkpoint_loaded": False,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
)
PY
rm -r "${RELOCATED}"
rm -r "${STAGE}/env"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
probe = json.loads((root / "validation" / "import_probe.json").read_text())
relocation = json.loads((root / "validation" / "relocation_probe.json").read_text())
if probe["compiled_mamba_imported"] or probe["checkpoint_loaded"] or probe["model_forward_executed"]:
    raise SystemExit("Caduceus runtime build executed model dependencies")
if (
    not relocation["prefix_relocated"]
    or relocation["compiled_mamba_imported"]
    or relocation["checkpoint_loaded"]
):
    raise SystemExit("Caduceus relocated runtime receipt differs")
freeze_tree(
    root,
    {
        "artifact_class": "caduceus_native_runtime",
        "upstream_source_artifacts_sha256": "4e1703b03ff28924904053e842b89079d510a80d2a9ad082de5369adaeea75c8",
        "python": "3.8.18",
        "pytorch": "2.2.0",
        "pytorch_cuda": "12.1",
        "mkl": "2023.1.0=h213fc3f_46344",
        "mamba_ssm": "1.2.0.post1+cu122torch2.2cxx11abiFALSE",
        "causal_conv1d": "1.2.0.post2+cu122torch2.2cxx11abiFALSE",
        "requirements_lock_sha256": sha256_file(root / "requirements.lock"),
        "runtime_payload_sha256": sha256_file(
            root / "payload" / "caduceus-native-env.tar.gz"
        ),
        "runtime_relocation_import_tested": True,
        "runtime_network_used_for_package_acquisition_only": True,
        "dependency_terms_review": "required_before_checkpoint_execution",
        "checkpoint_loaded": False,
    },
)
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Caduceus native runtime: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
