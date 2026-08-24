#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-acquire
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=executions/masld-corgi-acquire-%j.out
#SBATCH --error=executions/masld-corgi-acquire-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
CORGI_REV=d84c00082da3d04eb05169e290e53f1c7da2b9c5
REPRO_REV=2367c42cf36630b7cf561eccc2e781abc22b36be
ZENODO_RECORD=18630048
TARGET="${ROOT}/executions/upstream-acquisitions/corgi-d84c000-zenodo18630048-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite Corgi acquisition: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/downloads" "${STAGE}/source" "${STAGE}/validation"
cd "${ROOT}"

module purge
module load python/3.11.5-GCCcore-13.2.0

curl --fail --location --retry 5 --retry-all-errors --silent --show-error \
  "https://zenodo.org/api/records/${ZENODO_RECORD}" \
  --output "${STAGE}/downloads/zenodo-${ZENODO_RECORD}.json"

python - "${STAGE}/downloads/zenodo-${ZENODO_RECORD}.json" <<'PY'
from pathlib import Path
import json
import sys

record = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
if record.get("id") != 18630048:
    raise SystemExit("Zenodo Corgi record identity differs")
license_record = record.get("metadata", {}).get("license", {})
license_id = license_record.get("id") if isinstance(license_record, dict) else license_record
if str(license_id).lower() != "cc-by-4.0":
    raise SystemExit("Zenodo Corgi record license differs")
expected = {
    "corgi_model.pt": (2350451066, "md5:b10ddca8bb32eae3b67116242cfe1871"),
    "corgiplus_model.pt": (2384030158, "md5:f7dbe2ca301fe64814bebcd4fc3377bd"),
}
observed = {
    item["key"]: (item["size"], item["checksum"])
    for item in record.get("files", [])
    if item.get("key") in expected
}
if observed != expected:
    raise SystemExit("Zenodo Corgi checkpoint metadata differs")
PY

for filename in corgi_model.pt corgiplus_model.pt; do
  curl --fail --location --retry 5 --retry-all-errors --silent --show-error \
    "https://zenodo.org/records/${ZENODO_RECORD}/files/${filename}?download=1" \
    --output "${STAGE}/downloads/${filename}"
done

test "$(stat -c '%s' "${STAGE}/downloads/corgi_model.pt")" -eq 2350451066
test "$(stat -c '%s' "${STAGE}/downloads/corgiplus_model.pt")" -eq 2384030158
printf '%s  %s\n' b10ddca8bb32eae3b67116242cfe1871 \
  "${STAGE}/downloads/corgi_model.pt" | md5sum --check --strict
printf '%s  %s\n' f7dbe2ca301fe64814bebcd4fc3377bd \
  "${STAGE}/downloads/corgiplus_model.pt" | md5sum --check --strict

curl --fail --location --retry 5 --retry-all-errors --silent --show-error \
  "https://codeload.github.com/ekinda/corgi/tar.gz/${CORGI_REV}" \
  --output "${STAGE}/downloads/corgi-${CORGI_REV}.tar.gz"
curl --fail --location --retry 5 --retry-all-errors --silent --show-error \
  "https://codeload.github.com/ekinda/corgi-reproduction/tar.gz/${REPRO_REV}" \
  --output "${STAGE}/downloads/corgi-reproduction-${REPRO_REV}.tar.gz"

python - \
  "${STAGE}/downloads/corgi-${CORGI_REV}.tar.gz" \
  "${STAGE}/source" "corgi-${CORGI_REV}" \
  "${STAGE}/downloads/corgi-reproduction-${REPRO_REV}.tar.gz" \
  "${STAGE}/source" "corgi-reproduction-${REPRO_REV}" <<'PY'
from pathlib import Path, PurePosixPath
import sys
import tarfile

for offset in (1, 4):
    archive = Path(sys.argv[offset])
    destination = Path(sys.argv[offset + 1])
    expected_root = sys.argv[offset + 2]
    with tarfile.open(archive, "r:gz") as handle:
        members = handle.getmembers()
        if not members:
            raise SystemExit(f"empty source archive: {archive}")
        safe_members = []
        skipped_links = []
        for member in members:
            path = PurePosixPath(member.name)
            if path.is_absolute() or not path.parts or ".." in path.parts:
                raise SystemExit(f"unsafe source member: {member.name}")
            if path.parts[0] != expected_root:
                raise SystemExit(f"unexpected source archive root: {member.name}")
            if member.isdev() or member.islnk():
                raise SystemExit(f"device/hard-link source member rejected: {member.name}")
            if member.issym():
                expected_link = f"{expected_root}/corgi"
                expected_target = "/project/deeprna_data/corgi"
                if (
                    expected_root.startswith("corgi-reproduction-")
                    and member.name == expected_link
                    and member.linkname == expected_target
                ):
                    skipped_links.append({"path": member.name, "target": member.linkname})
                    continue
                raise SystemExit(f"unexpected source symlink rejected: {member.name}")
            safe_members.append(member)
        if expected_root.startswith("corgi-reproduction-") and skipped_links != [{
            "path": f"{expected_root}/corgi",
            "target": "/project/deeprna_data/corgi",
        }]:
            raise SystemExit("expected nonportable Corgi reproduction symlink not found")
        if not expected_root.startswith("corgi-reproduction-") and skipped_links:
            raise SystemExit("unexpected skipped source links")
        handle.extractall(destination, members=safe_members)
PY

python - \
  "${STAGE}/source/corgi-${CORGI_REV}" \
  "${STAGE}/source/corgi-reproduction-${REPRO_REV}" \
  "${STAGE}/validation/source-contract.json" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

corgi = Path(sys.argv[1])
reproduction = Path(sys.argv[2])
expected = {
    corgi / "LICENSE": (9707, "feb79b816ab6d0cfef7c005b934f1383c0c1af9b6522382a2e156b82f87e388a"),
    corgi / "pyproject.toml": (725, "fd298c3d18d273d8b41fa6df57a4f1d52cba4f7b2f68e35a57537090f66ef515"),
    corgi / "corgi/config.py": (2256, "a6ae55409602e0de35a8d2cf757eb5ebac93b602e31e00f0aa36a4fbd63ac4b1"),
    corgi / "corgi/constants.py": (674, "4ee89ba8fef8c72f931dabd9b25545e95922b9c1a9d1da8235f951b99b39ad7c"),
    corgi / "corgi/model.py": (12228, "767067e4659ca196f307d2133104365fb7e6ea9d2f151dd61d3c2a8a7e0737e6"),
    corgi / "corgi/predict.py": (17976, "0d8da0bedbf60a493188bfd9109cc56624292ea139c14f3ec20df5cf8573567b"),
    corgi / "data/tf_reference.npy": (23256, "95041201bacaa6a77f0cf8c6db29a628b47249edcc0ee69579d04d49b6b2381d"),
    corgi / "corgi/config_corgiplus.py": (2410, "9d996d85a8e1ab8b3b0a21277581a1ddb5efec65f9a28c841215986a1f9c13b8"),
    corgi / "corgi/trainer_corgiplus.py": (25752, "4ccb5db5d2be4933aa0cb270217fe251e9eb71e9635ed7cfcacb879c27926af5"),
    reproduction / "code/corgiplus_delta_bench_tta.py": (24502, "a6ec48da151f6da912d4b6d019a3e0d021322867c91ec7a5f6a920ccb9212276"),
}
records = []
for path, (expected_size, expected_hash) in expected.items():
    if not path.is_file() or path.is_symlink():
        raise SystemExit(f"missing or unsafe pinned source: {path}")
    observed_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    if path.stat().st_size != expected_size or observed_hash != expected_hash:
        raise SystemExit(f"pinned Corgi source differs: {path}")
    records.append({
        "identity": (
            f"ekinda/corgi@d84c000:{path.relative_to(corgi).as_posix()}"
            if path.is_relative_to(corgi)
            else "ekinda/corgi-reproduction@2367c42:"
            f"{path.relative_to(reproduction).as_posix()}"
        ),
        "size_bytes": expected_size,
        "sha256": expected_hash,
    })
license_text = (corgi / "LICENSE").read_text(encoding="utf-8")
pyproject_text = (corgi / "pyproject.toml").read_text(encoding="utf-8")
predict_text = (corgi / "corgi/predict.py").read_text(encoding="utf-8")
if "Apache License" not in license_text or 'license = {text = "MIT"}' not in pyproject_text:
    raise SystemExit("Corgi source license declarations differ")
if "raise NotImplementedError" not in predict_text or "weights_only=True" in predict_text:
    raise SystemExit("Corgi released loader audit assumptions differ")
Path(sys.argv[3]).write_text(json.dumps({
    "status": "pass",
    "regular_high_level_api_present": "class corgi_pretrained" in predict_text,
    "corgi_plus_placeholder_present": "class corgiplus_pretrained" in predict_text,
    "released_loader_weights_only": False,
    "code_license_conflict_present": True,
    "skipped_nonportable_reproduction_symlink": {
        "path": "corgi",
        "target": "/project/deeprna_data/corgi",
    },
    "records": records,
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

python - "${STAGE}/downloads" "${STAGE}/validation/checkpoints.json" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

root = Path(sys.argv[1])
expected = {
    "corgi_model.pt": (2350451066, "b10ddca8bb32eae3b67116242cfe1871"),
    "corgiplus_model.pt": (2384030158, "f7dbe2ca301fe64814bebcd4fc3377bd"),
}
records = []
for name, (size, md5) in expected.items():
    path = root / name
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            sha256.update(chunk)
    records.append({
        "filename": name,
        "size_bytes": size,
        "md5": md5,
        "sha256": sha256.hexdigest(),
        "deserialized": False,
    })
Path(sys.argv[2]).write_text(json.dumps({
    "schema_version": "masld-bench-corgi-checkpoint-acquisition-v1",
    "zenodo_record": 18630048,
    "record_license": "CC-BY-4.0",
    "checkpoints": records,
}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
PY

cat > "${STAGE}/SOURCE.tsv" <<EOF
field\tvalue
model_family\tcorgi
code_repository\tekinda/corgi
code_revision\t${CORGI_REV}
reproduction_repository\tekinda/corgi-reproduction
reproduction_revision\t${REPRO_REV}
zenodo_record\t${ZENODO_RECORD}
checkpoint_deserialization\tfalse
EOF

sha256sum \
  "${ROOT}/config/artifacts/models/corgi/checkpoints.json" \
  "${ROOT}/config/artifacts/models/corgi_plus/checkpoints.json" \
  "${ROOT}/config/artifacts/models/corgi/release.json" \
  "${ROOT}/slurm/corgi_acquire_official_release.sbatch" \
  > "${STAGE}/validation/input-locks.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

freeze_tree(
    Path(sys.argv[1]),
    metadata={
        "artifact_class": "corgi_official_release_acquisition",
        "model_ids": ["corgi_regular", "corgi_plus"],
        "code_revision": "d84c00082da3d04eb05169e290e53f1c7da2b9c5",
        "reproduction_revision": "2367c42cf36630b7cf561eccc2e781abc22b36be",
        "zenodo_record": 18630048,
        "checkpoint_deserialized": False,
        "evaluator_outcomes_exposed": False,
    },
)
PY

mv "${STAGE}" "${TARGET}"
printf 'Corgi official release acquisition: %s\n' "${TARGET}"
