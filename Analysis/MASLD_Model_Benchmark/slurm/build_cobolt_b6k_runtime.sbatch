#!/usr/bin/env bash
#SBATCH --job-name=model-env-304
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# Build the Cobolt B6K runtime: a fresh micromamba environment whose ONLY
# scientific delta from the validated L40S runtime (.agent/scanpy_env,
# environment pin gpu_rna_atac_torch_smoke--c2308fc3b7524081) is torch, which
# moves 2.3.1+cu121 -> 2.8.0+cu128 so it can execute on sm_120.
#
# Python stays 3.10 and every package named in the L40S selection record's critical_versions
# is pinned to its exact L40S value, so the re-validation isolates the torch
# change.  .agent/scanpy_env is READ for its pin values and never modified.
#
# The environment is created directly at its final path.  A conda prefix embeds
# absolute paths in script shebangs, so a staging directory followed by a rename
# would corrupt it.  Completion is marked by removing the INCOMPLETE sentinel.
#
# freeze_tree() is NOT applied to this tree: inventory_tree() rejects symlinks and
# a conda prefix is full of them.  This follows the existing environment-tree
# convention (see scprint2-native-torch2.8.0-v2-21073597), where SHA256SUMS and
# ARTIFACTS.json cover the metadata files while the environment itself is pinned
# by pip_freeze_sha256 in the separate environment pin output file.
#
# Installs packages at BUILD time only.  Reads no biological data, no development
# outcome, no held-back outcome, and calculates no benchmark metric.  Writes only
# under executions/, which the fixed source record tree walker excludes.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
MICROMAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
L40S_LOCK="${ROOT}/executions/environment-locks/gpu_rna_atac_torch_smoke--c2308fc3b7524081"
L40S_LOCK_ARTIFACTS_SHA256=f67836fe3618a61289bd4efb5d64367c2ecf310ea10a441e62f660aa464a720c
L40S_LOCK_JSON_SHA256=c2308fc3b7524081c7a49ff1d840213306bcc20dff003010da2b898e2bdaa461
RUNTIME_ID=gpu_rna_atac_torch_b6k
TARGET="${ROOT}/executions/environments/cobolt-native-torch2.8.0-cu128-${SLURM_JOB_ID}"
FAILED="${TARGET}.failed"
ENV_PREFIX="${TARGET}/env"

preserve_failed_tree() {
  if [[ -e "${TARGET}/INCOMPLETE" && ! -e "${FAILED}" ]]; then mv "${TARGET}" "${FAILED}"; fi
}
trap preserve_failed_tree EXIT

if [[ -e "${TARGET}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite runtime tree: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${TARGET}"
: > "${TARGET}/INCOMPLETE"

cd "${ROOT}"

# The L40S selection record is the source of every pin below; prove it has not drifted.
sha256sum --check --strict <<EOF
${L40S_LOCK_ARTIFACTS_SHA256}  ${L40S_LOCK}/ARTIFACTS.json
${L40S_LOCK_JSON_SHA256}  ${L40S_LOCK}/environment.lock.json
EOF

"${MICROMAMBA}" --version > "${TARGET}/micromamba_version.txt"

python3 - "${L40S_LOCK}/environment.lock.json" "${TARGET}/l40s_reference_pins.json" <<'PYA'
import json, sys
lock = json.loads(open(sys.argv[1]).read())
expected = {
    "anndata": "0.11.4", "h5py": "3.15.1", "numpy": "2.2.6", "pandas": "2.3.3",
    "scanpy": "1.11.5", "scikit-learn": "1.7.2", "scipy": "1.15.3",
}
if lock["critical_versions"] != expected:
    raise SystemExit(f"L40S critical_versions drifted: {lock['critical_versions']}")
if not lock["python_version"].startswith("3.10."):
    raise SystemExit(f"L40S python is not 3.10.x: {lock['python_version']}")
json.dump(
    {
        "critical_versions": expected,
        "l40s_python_version": lock["python_version"],
        "torch_after": "2.8.0+cu128",
        "torch_before": "2.3.1+cu121",
    },
    open(sys.argv[2], "w"), indent=2, sort_keys=True,
)
PYA

export CONDA_PKGS_DIRS="${TMPDIR:-/tmp}/cobolt-b6k-pkgs-${SLURM_JOB_ID}"

# conda layer: exactly the L40S critical_versions, except scipy.  conda-forge has
# no scipy 1.15.3 for py310; the L40S env carries conda scipy 1.15.2 shadowed by
# pip scipy 1.15.3, and that two-layer shape is reproduced below.
"${MICROMAMBA}" create --yes --no-rc --prefix "${ENV_PREFIX}" -c conda-forge \
  'python=3.10.19' pip \
  'numpy=2.2.6' 'h5py=3.15.1' 'pandas=2.3.3' 'anndata=0.11.4' \
  'scanpy=1.11.5' 'scikit-learn=1.7.2' 'scipy=1.15.2' \
  'matplotlib-base=3.9.3' 'seaborn-base=0.13.2'

PY_BIN="${ENV_PREFIX}/bin/python"

# pip layer, one deterministic index per call so there is no cross-index race.
"${PY_BIN}" -m pip install --no-cache-dir --no-input \
  --index-url https://download.pytorch.org/whl/cu128 \
  'torch==2.8.0+cu128'
"${PY_BIN}" -m pip install --no-cache-dir --no-input \
  --index-url https://pypi.org/simple \
  'scipy==1.15.3'

"${PY_BIN}" -m pip freeze --all > "${TARGET}/pip_freeze.txt"
"${MICROMAMBA}" list --explicit --sha256 --no-rc --prefix "${ENV_PREFIX}" \
  > "${TARGET}/conda_explicit.txt"

"${PY_BIN}" - <<'PYB' > "${TARGET}/runtime_probe.json"
import importlib.metadata as md
import json, os, platform, subprocess, sys

import torch
import matplotlib, seaborn  # upstream coboltmodel.py imports these at module level

expected = {
    "anndata": "0.11.4", "h5py": "3.15.1", "numpy": "2.2.6", "pandas": "2.3.3",
    "scanpy": "1.11.5", "scikit-learn": "1.7.2", "scipy": "1.15.3",
}
observed = {name: md.version(name) for name in expected}
if observed != expected:
    raise SystemExit(f"critical versions do not match the L40S runtime: {observed}")
if torch.__version__ != "2.8.0+cu128":
    raise SystemExit(f"torch is {torch.__version__}")
if not platform.python_version().startswith("3.10."):
    raise SystemExit(f"python is {platform.python_version()}")
if torch.version.cuda != "12.8":
    raise SystemExit(f"torch cuda runtime is {torch.version.cuda}")
# get_arch_list() needs an initialised driver and returns [] on this GPU-less io
# node.  sm_120 SASS is asserted on real hardware by the B6K validation job.
arch = list(torch.cuda.get_arch_list())

# base._validate_environment shells out to this; the scprint2 cu128 venv fails here.
pip = subprocess.run(
    [sys.executable, "-m", "pip", "freeze", "--all"],
    check=False, capture_output=True, text=True,
    env={**os.environ, "PYTHONNOUSERSITE": "1"},
)
if pip.returncode != 0:
    raise SystemExit("pip freeze --all must succeed for base._validate_environment")

json.dump({
    "arch_list_on_build_node": arch,
    "critical_versions": observed,
    "cuda_available_on_build_node": torch.cuda.is_available(),
    "cuda_runtime": torch.version.cuda,
    "environment_prefix": os.path.realpath(sys.prefix),
    "matplotlib": matplotlib.__version__,
    "pip_freeze_available": True,
    "platform": platform.platform(),
    "python": sys.version,
    "python_executable": os.path.realpath(sys.executable),
    "python_version": platform.python_version(),
    "runtime_id": "gpu_rna_atac_torch_b6k",
    "seaborn": seaborn.__version__,
    "torch": torch.__version__,
    "torch_file": os.path.realpath(torch.__file__),
}, sys.stdout, indent=2, sort_keys=True)
print()
PYB

# Separate environment pin output file: the record base._validate_environment reads.
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="${ROOT}/src" python "${ROOT}/scripts/freeze_python_environment.py" \
  --micromamba "${MICROMAMBA}" \
  --prefix "${ENV_PREFIX}" \
  --runtime-id "${RUNTIME_ID}" \
  --output-root "${ROOT}/executions/environment-locks" \
  > "${TARGET}/environment_lock_record.json"
cat "${TARGET}/environment_lock_record.json"

# Metadata manifest.  Deliberately excludes env/, per the note in the header.
python3 - "${TARGET}" "${RUNTIME_ID}" <<'PYC'
import hashlib, json, sys
from pathlib import Path

target, runtime_id = Path(sys.argv[1]), sys.argv[2]

def digest(name):
    return hashlib.sha256((target / name).read_bytes()).hexdigest()

record = json.loads((target / "environment_lock_record.json").read_text())
payload = {
    "conda_explicit_sha256": digest("conda_explicit.txt"),
    "environment_lock_directory": record["lock_directory"],
    "environment_lock_sha256": record["environment_artifact"]["sha256"],
    "environment_prefix": (target / "env").as_posix(),
    "environment_tree_excludes_env_directory": True,
    "model_id": "cobolt",
    "package_installation_performed_at_build_only": True,
    "pip_freeze_sha256": digest("pip_freeze.txt"),
    "runtime_id": runtime_id,
    "runtime_probe_sha256": digest("runtime_probe.json"),
    "schema_version": "masld-bench-cobolt-b6k-runtime-v1",
    "status": "pass",
    "torch_required": "2.8.0+cu128",
    "upstream_delta_from_l40s_runtime": "torch only",
}
(target / "ARTIFACTS.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(json.dumps(payload, indent=2, sort_keys=True))
PYC

(
  cd "${TARGET}"
  sha256sum ARTIFACTS.json conda_explicit.txt environment_lock_record.json \
    l40s_reference_pins.json micromamba_version.txt pip_freeze.txt \
    runtime_probe.json > SHA256SUMS
)

rm -f "${TARGET}/INCOMPLETE"
trap - EXIT
printf 'Cobolt B6K runtime built: %s\n' "${TARGET}"
