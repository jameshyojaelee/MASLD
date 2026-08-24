#!/usr/bin/env bash
#SBATCH --job-name=masld-evo-vtx-admit
#SBATCH --partition=io
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
EVO_REV=53f195997257c56c00e5ef8d33a54f5baad143a6
TARGET="${ROOT}/executions/evo-vtx-runtime-admission-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
    printf 'Preserved failed Evo/Vortex runtime admission: %s\n' "${FAILED}" >&2
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Evo/Vortex admission: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/downloads" "${STAGE}/metadata" \
  "${STAGE}/inventories" "${STAGE}/validation"
curl_flags=(--fail --location --retry 6 --retry-all-errors --connect-timeout 30 \
  --proto '=https' --tlsv1.2)
curl "${curl_flags[@]}" --output "${STAGE}/metadata/vtx-1.1.0.pypi.json" \
  https://pypi.org/pypi/vtx/1.1.0/json
curl "${curl_flags[@]}" --output "${STAGE}/downloads/evo2-source.tar.gz" \
  "https://codeload.github.com/ArcInstitute/evo2/tar.gz/${EVO_REV}"

module load python/3.11.5-GCCcore-13.2.0
mapfile -t vtx_urls < <(python - "${STAGE}/metadata/vtx-1.1.0.pypi.json" <<'PY'
import json
import sys
data = json.load(open(sys.argv[1]))
expected = {
    "vtx-1.1.0-py3-none-any.whl": "0ff9f1db2f9e81e288150b60fd4fe4832b8b992ac2c6c947271b2036ffeb8299",
    "vtx-1.1.0.tar.gz": "0d54516d8a76aeec6839539b17273fa763349041266b6be7a1ae85c59f6d2404",
}
observed = {
    row["filename"]: row
    for row in data["urls"]
    if row["filename"] in expected
}
if set(observed) != set(expected) or any(
    observed[name]["digests"]["sha256"] != digest
    for name, digest in expected.items()
):
    raise SystemExit("Vortex PyPI release identity differs")
for name in sorted(expected):
    print(name + "\t" + observed[name]["url"])
PY
)
if [[ ${#vtx_urls[@]} != 2 ]]; then
  printf 'Vortex URL census differs\n' >&2
  exit 1
fi
for binding in "${vtx_urls[@]}"; do
  IFS=$'\t' read -r name url <<< "${binding}"
  curl "${curl_flags[@]}" --output "${STAGE}/downloads/${name}" "${url}"
done
printf '%s  %s\n%s  %s\n' \
  0ff9f1db2f9e81e288150b60fd4fe4832b8b992ac2c6c947271b2036ffeb8299 "${STAGE}/downloads/vtx-1.1.0-py3-none-any.whl" \
  0d54516d8a76aeec6839539b17273fa763349041266b6be7a1ae85c59f6d2404 "${STAGE}/downloads/vtx-1.1.0.tar.gz" \
  | sha256sum --check --strict

python "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  --archive "${STAGE}/downloads/vtx-1.1.0.tar.gz" \
  --output "${STAGE}/inventories/vtx-source.inventory.json" \
  > "${STAGE}/validation/vtx-source.inventory.stdout.txt"
python "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  --archive "${STAGE}/downloads/evo2-source.tar.gz" \
  --output "${STAGE}/inventories/evo2-source.inventory.json" \
  > "${STAGE}/validation/evo2-source.inventory.stdout.txt"
python - "${STAGE}" "${EVO_REV}" <<'PY'
from hashlib import sha256
import io
import json
from pathlib import Path, PurePosixPath
import stat
import sys
import tarfile
import zipfile

root = Path(sys.argv[1])
revision = sys.argv[2]
with zipfile.ZipFile(root / "downloads/vtx-1.1.0-py3-none-any.whl") as archive:
    members = []
    seen = set()
    total = 0
    for info in archive.infolist():
        path = PurePosixPath(info.filename)
        mode = info.external_attr >> 16
        if (
            not info.filename
            or "\\" in info.filename
            or path.is_absolute()
            or ".." in path.parts
            or info.filename in seen
            or stat.S_ISLNK(mode)
            or info.flag_bits & 1
        ):
            raise SystemExit("Vortex wheel contains an unsafe member")
        seen.add(info.filename)
        if info.is_dir():
            members.append({"path": info.filename, "type": "directory", "size_bytes": 0})
            continue
        total += info.file_size
        if total > 1_000_000_000:
            raise SystemExit("Vortex wheel uncompressed size differs")
        digest = sha256()
        observed = 0
        with archive.open(info) as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
                observed += len(block)
        if observed != info.file_size:
            raise SystemExit("Vortex wheel member size differs")
        members.append({
            "path": info.filename,
            "type": "regular_file",
            "size_bytes": observed,
            "sha256": digest.hexdigest(),
        })
    (root / "inventories/vtx-wheel.inventory.json").write_text(json.dumps({
        "schema_version": "masld-bench-safe-wheel-inventory-v1",
        "status": "pass",
        "archive_extracted": False,
        "member_count": len(members),
        "uncompressed_regular_bytes": total,
        "members": members,
    }, sort_keys=True, indent=2) + "\n")
    metadata_names = [name for name in archive.namelist() if name.endswith(".dist-info/METADATA")]
    if len(metadata_names) != 1:
        raise SystemExit("Vortex wheel metadata census differs")
    metadata = archive.read(metadata_names[0]).decode("utf-8")
requirements = sorted(
    line.split(":", 1)[1].strip()
    for line in metadata.splitlines()
    if line.startswith("Requires-Dist:")
)
if "Version: 1.1.0" not in metadata or "Name: vtx" not in metadata:
    raise SystemExit("Vortex wheel identity differs")

expected_sources = {
    "LICENSE": "5bb5812fc2bfb2d777fe5621767172f8ef30be62ad719c87e0f10832336a99e0",
    "NOTICE": "70672511c8a6ac6eadc2f626446a47c8a72c577326262f99244571be72f51036",
    "README.md": "58787c8ef5cb4fba4c04322a4ceb9f174e2233ec22d4193622fb6bc67d651d89",
    "pyproject.toml": "7cc6a516e1da69df2f87f6487f80bdc69a50b581c4a674a3313b4c9c23b5cddf",
    "evo2/configs/evo2-7b-1m.yml": "9637f84ffd796f42cce3f926832b3999ff34a82686beba04505327162b288b1a",
    "evo2/scoring.py": "17e6985cef3a62c2a488133ab5a65050dcf0abeec2f55d1c0b068dfe37c44de6",
    "evo2/models.py": "0b757fcc42b37a52b0c3f508e876a52356331392ae8aa16a923576265043b78c",
}
with tarfile.open(root / "downloads/evo2-source.tar.gz", "r:gz") as archive:
    prefix = f"evo2-{revision}/"
    for relative, expected in expected_sources.items():
        member = archive.getmember(prefix + relative)
        stream = archive.extractfile(member)
        if stream is None or sha256(stream.read()).hexdigest() != expected:
            raise SystemExit(f"Evo source identity differs: {relative}")

receipt = {
    "schema_version": "masld-bench-evo-vtx-runtime-admission-v1",
    "status": "pass",
    "evo2_source_revision": revision,
    "vtx_version": "1.1.0",
    "vtx_wheel_sha256": "0ff9f1db2f9e81e288150b60fd4fe4832b8b992ac2c6c947271b2036ffeb8299",
    "vtx_source_sha256": "0d54516d8a76aeec6839539b17273fa763349041266b6be7a1ae85c59f6d2404",
    "vtx_requires_dist": requirements,
    "archives_extracted": False,
    "custom_model_code_imported": False,
    "model_instantiated": False,
    "checkpoint_loaded": False,
    "observed_outcomes_loaded": False,
    "sealed_outcomes_loaded": False,
}
(root / "validation/admission_receipt.json").write_text(
    json.dumps(receipt, sort_keys=True, indent=2) + "\n"
)
PY
python -m pip freeze > "${STAGE}/validation/effective-pip-freeze.txt"
sha256sum \
  "${ROOT}/scripts/epibert_safe_archive_inventory.py" \
  "${ROOT}/slurm/evo_vtx_runtime_admission.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "validation/admission_receipt.json").read_text())
if (
    receipt["status"] != "pass"
    or receipt["vtx_version"] != "1.1.0"
    or receipt["archives_extracted"]
    or receipt["custom_model_code_imported"]
    or receipt["model_instantiated"]
    or receipt["checkpoint_loaded"]
    or receipt["observed_outcomes_loaded"]
    or receipt["sealed_outcomes_loaded"]
):
    raise SystemExit("Evo/Vortex runtime admission receipt differs")
freeze_tree(root, {
    "artifact_class": "evo_vtx_runtime_source_admission",
    "evo2_source_revision": receipt["evo2_source_revision"],
    "vtx_version": receipt["vtx_version"],
    "vtx_wheel_sha256": receipt["vtx_wheel_sha256"],
    "vtx_source_sha256": receipt["vtx_source_sha256"],
    "terms": "Apache-2.0_with_NOTICE_dependency_audit_pending",
    "checkpoint_loaded": False,
    "observed_outcomes_loaded": False,
    "sealed_outcomes_loaded": False,
})
PY
module purge

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Evo/Vortex runtime admission: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
