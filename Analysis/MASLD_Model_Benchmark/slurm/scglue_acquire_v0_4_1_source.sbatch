#!/usr/bin/env bash
#SBATCH --job-name=masld-scglue041-source
#SBATCH --partition=io
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
REVISION=091bdbbdb3dd6d652ed4677200b034ee3c85dae1
URL="https://codeload.github.com/gao-lab/GLUE/tar.gz/${REVISION}"
TARGET="${ROOT}/executions/upstream-acquisitions/scglue-v0.4.1-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

preserve_failed_stage() {
  local failed="${TARGET}.failed"
  if [[ -d "${STAGE}" ]]; then
    if [[ -e "${failed}" ]]; then
      printf 'Cannot preserve failed scGLUE source stage: %s\n' "${failed}" >&2
      return
    fi
    mv "${STAGE}" "${failed}"
    printf 'Preserved failed scGLUE source stage: %s\n' "${failed}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${TARGET}.failed" ]]; then
  printf 'Refusing to overwrite scGLUE source acquisition: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/download" "${STAGE}/source" "${STAGE}/validation"
curl --fail --location --retry 4 --retry-all-errors \
  --output "${STAGE}/download/scglue-${REVISION}.tar.gz" \
  "${URL}"

module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}" "${REVISION}" <<'PY'
from pathlib import Path, PurePosixPath
from hashlib import sha256
import json
import sys
import tarfile

root = Path(sys.argv[1]).resolve(strict=True)
revision = sys.argv[2]
archive = root / f"download/scglue-{revision}.tar.gz"
expected = {
    "LICENSE": ("f66904d80415738a11a60a542ebe29688993481f415cbc82f39a1cb7ed6c82fe", 1064),
    "pyproject.toml": ("e51acb7bc8975176e8f7cba39d6b193118c5dc1f6674abfced9c05b4c9ed2604", 2048),
    "scglue/models/scglue.py": ("00343c0aa394f23858eb1ba6d526bd3c61ec6d5a33c2c0e2e99d645c442ace84", 48463),
    "scglue/models/base.py": ("5b507022f5bd48148ef59647038bf1d8169a0eb7cdb908670a505afe14570b2b", 11827),
    "scglue/genomics.py": ("285b4e4bdb32223575241aec3568faf83c1b216498b354439b7dea4ee7eec3d4", 39165),
}
with tarfile.open(archive, "r:gz") as handle:
    members = handle.getmembers()
    if not members or len(members) > 10000:
        raise SystemExit("scGLUE archive member count differs")
    roots = set()
    total = 0
    safe_members = []
    excluded_symlinks = []
    for member in members:
        path = PurePosixPath(member.name)
        if path.is_absolute() or ".." in path.parts:
            raise SystemExit(f"unsafe scGLUE archive member: {member.name!r}")
        roots.add(path.parts[0])
        if member.issym():
            target = PurePosixPath(member.linkname)
            if target.is_absolute():
                raise SystemExit(f"unsafe scGLUE symlink target: {member.linkname!r}")
            normalized = []
            for part in (path.parent / target).parts:
                if part in {"", "."}:
                    continue
                if part == "..":
                    if not normalized:
                        raise SystemExit(f"escaping scGLUE symlink target: {member.linkname!r}")
                    normalized.pop()
                else:
                    normalized.append(part)
            if not normalized or normalized[0] != path.parts[0]:
                raise SystemExit(f"escaping scGLUE symlink target: {member.linkname!r}")
            excluded_symlinks.append({
                "path": member.name,
                "target": member.linkname,
                "normalized_target": "/".join(normalized),
            })
            continue
        if not (member.isdir() or member.isfile()):
            raise SystemExit(f"unsafe scGLUE archive member: {member.name!r}")
        safe_members.append(member)
        total += member.size
    if len(roots) != 1 or total > 500_000_000:
        raise SystemExit("scGLUE archive root or expansion differs")
    source_root = next(iter(roots))
    handle.extractall(root / "source", members=safe_members, filter="data")

resolved = root / "source" / source_root
records = {}
for relative, (digest_expected, size_expected) in expected.items():
    path = resolved / relative
    digest = sha256(path.read_bytes()).hexdigest()
    if digest != digest_expected or path.stat().st_size != size_expected:
        raise SystemExit(f"scGLUE source identity differs: {relative}")
    records[relative] = {"sha256": digest, "size_bytes": path.stat().st_size}
(root / "validation/source-admission.json").write_text(
    json.dumps({
        "schema_version": "masld-bench-scglue-source-admission-v1",
        "status": "pass",
        "repository": "gao-lab/GLUE",
        "revision": revision,
        "tag": "v0.4.1",
        "source_root": source_root,
        "archive_sha256": sha256(archive.read_bytes()).hexdigest(),
        "archive_size_bytes": archive.stat().st_size,
        "expanded_size_bytes": total,
        "excluded_safe_relative_symlinks": excluded_symlinks,
        "critical_source_records": records,
        "code_license": "MIT",
        "checkpoint_downloaded": False,
        "serialized_object_loaded": False,
        "project_data_read": False,
        "outcomes_read": False,
    }, sort_keys=True, indent=2) + "\n",
    encoding="utf-8",
)
PY
sha256sum \
  "${ROOT}/slurm/scglue_acquire_v0_4_1_source.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" "${REVISION}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/source-admission.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["revision"] != sys.argv[2]
    or receipt["code_license"] != "MIT"
    or receipt["checkpoint_downloaded"]
    or receipt["serialized_object_loaded"]
    or receipt["project_data_read"]
    or receipt["outcomes_read"]
):
    raise SystemExit("scGLUE source admission differs")
freeze_tree(
    root,
    {
        "artifact_class": "scglue_v0_4_1_source",
        "repository": "gao-lab/GLUE",
        "revision": sys.argv[2],
        "code_license": "MIT",
        "project_data_read": False,
        "outcomes_read": False,
        "runtime_probe_pending": True,
        "profile_task_status": "latent_only",
        "status": "passed",
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'scGLUE v0.4.1 source: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
