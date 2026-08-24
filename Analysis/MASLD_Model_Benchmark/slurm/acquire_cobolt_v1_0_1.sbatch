#!/usr/bin/env bash
#SBATCH --job-name=masld-cobolt-acquire
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/tests/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
OUT="${ROOT}/executions/upstream-acquisitions/cobolt-v1.0.1-${SLURM_JOB_ID}"
ARCHIVE="${OUT}/epurdom-cobolt-v1.0.1.zip"
URL='https://zenodo.org/api/records/7039226/files/epurdom/cobolt-v1.0.1.zip/content'
EXPECTED_SHA256=ee72da3d1869d31d800441c021543e8ad4af3fcc99ca0cf630197801aba632d0
EXPECTED_MD5=7341b3a8d1f8a05a082029499c6c0af0

mkdir -p "${OUT}"
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0

curl --fail --location --proto '=https' --tlsv1.2 \
  --retry 3 --retry-all-errors --connect-timeout 30 \
  --output "${ARCHIVE}" "${URL}"
printf '%s  %s\n' "${EXPECTED_SHA256}" "${ARCHIVE}" | sha256sum -c -
printf '%s  %s\n' "${EXPECTED_MD5}" "${ARCHIVE}" | md5sum -c -

python - "${ARCHIVE}" "${OUT}" <<'PY'
from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import stat
import sys
import zipfile

archive = Path(sys.argv[1]).resolve(strict=True)
out = Path(sys.argv[2]).resolve(strict=True)
source = out / "source"
expected = {
    "LICENSE": "3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986",
    "README.md": "64d15e5afb832c9d00ebe91be104807c9be010403a361458a82e25aef950c340",
    "setup.py": "71ae102c492515c2e047323781243b12837d32b4dd5811ade6669d7c37cf0204",
    "cobolt/model/cobolt.py": "6481b6a89d3edde3c9b51d1d3ffb36a58e8762818cc30d585601a457b2514001",
    "cobolt/model/coboltmodel.py": "a4f5c3a45ca4069e12ffa05d3efbfcd8a0a636da7515d7d4e9174a914f3f2311",
    "cobolt/utils/data.py": "f2d7937753fd3c92aecd34ebf56913b957f862caf167f990eac3a71c413e7aa4",
    "cobolt/utils/dataset.py": "8fa2b8021f0cc5980ccc3a4625a9ea43e444a105a4797c5f3cc5a94c97a6ad88",
}

with zipfile.ZipFile(archive) as handle:
    infos = handle.infolist()
    if not infos or len(infos) > 10_000:
        raise RuntimeError("Cobolt archive entry count is invalid")
    total_size = 0
    inventory = []
    for info in infos:
        path = PurePosixPath(info.filename)
        mode = (info.external_attr >> 16) & 0xFFFF
        if (
            path.is_absolute()
            or ".." in path.parts
            or not path.parts
            or info.flag_bits & 0x1
            or stat.S_ISLNK(mode)
            or info.file_size > 50_000_000
        ):
            raise RuntimeError(f"unsafe Cobolt ZIP entry: {info.filename}")
        total_size += info.file_size
        inventory.append(
            {
                "path": info.filename,
                "compressed_size": info.compress_size,
                "size_bytes": info.file_size,
            }
        )
    if total_size > 500_000_000:
        raise RuntimeError("Cobolt archive expands beyond the admission limit")
    roots = {PurePosixPath(info.filename).parts[0] for info in infos}
    if len(roots) != 1:
        raise RuntimeError("Cobolt archive does not have one source root")
    source_root = source / next(iter(roots))
    source.mkdir()
    handle.extractall(source)

matches = {}
for suffix, expected_sha in expected.items():
    path = source_root / suffix
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"Cobolt source identity is missing: {suffix}")
    observed = sha256(path.read_bytes()).hexdigest()
    if observed != expected_sha:
        raise RuntimeError(f"Cobolt source hash differs: {suffix}")
    matches[suffix] = {
        "path": path.relative_to(out).as_posix(),
        "sha256": observed,
        "size_bytes": path.stat().st_size,
    }

payload = {
    "schema_version": "masld-bench-safe-upstream-source-acquisition-v1",
    "repository": "epurdom/cobolt",
    "repository_revision": "cf5a448c6539025346a6393c215ba329cb3a0183",
    "release_tag": "v1.0.1",
    "zenodo_record": "7039226",
    "archive_entry_count": len(inventory),
    "archive_expanded_size_bytes": total_size,
    "unsafe_entries_detected": False,
    "source_hashes_verified": matches,
    "inventory": inventory,
}
(out / "safe_inventory.json").write_text(
    json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n",
    encoding="utf-8",
)
PY

PYTHONPATH=src python - "${OUT}" "${EXPECTED_SHA256}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
archive = root / "epurdom-cobolt-v1.0.1.zip"
freeze_tree(
    root,
    {
        "artifact_class": "upstream_source_acquisition",
        "archive_sha256": sha256_file(archive),
        "expected_archive_sha256": sys.argv[2],
        "model_id": "cobolt",
        "release_tag": "v1.0.1",
        "source_revision": "cf5a448c6539025346a6393c215ba329cb3a0183",
    },
)
print(root)
print(sha256_file(root / "ARTIFACTS.json"))
PY
