#!/usr/bin/env bash
#SBATCH --job-name=source-fetch-256398
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
BASE_URL=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE256nnn/GSE256398
ARCHIVE_URL=${BASE_URL}/suppl/GSE256398_RAW.tar
FILELIST_URL=${BASE_URL}/suppl/filelist.txt
SOFT_URL=${BASE_URL}/soft/GSE256398_family.soft.gz
TARGET=${ROOT}/executions/gse256398-source-h5-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite GSE256398 acquisition: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/archive" "${STAGE}/h5" "${STAGE}/metadata" "${STAGE}/validation"

curl_flags=(--fail --location --retry 8 --retry-all-errors --connect-timeout 30)
curl "${curl_flags[@]}" --output "${STAGE}/metadata/filelist.txt" "${FILELIST_URL}"
curl "${curl_flags[@]}" --output "${STAGE}/metadata/GSE256398_family.soft.gz" "${SOFT_URL}"
test "$(stat -c '%s' "${STAGE}/metadata/filelist.txt")" -eq 2764
test "$(stat -c '%s' "${STAGE}/metadata/GSE256398_family.soft.gz")" -eq 4886
printf '%s  %s\n' \
  ff28ae965f129f9014043f97fc9811007a89692ad89ab5e2f46fff3903ecd674 "${STAGE}/metadata/filelist.txt" \
  5d5b354b283d7a92693ce604037b3e0eac5b01d6da0a5b51c9bdc93be8655add "${STAGE}/metadata/GSE256398_family.soft.gz" \
  | sha256sum --check --strict
gzip -t "${STAGE}/metadata/GSE256398_family.soft.gz"

ARCHIVE=${STAGE}/archive/GSE256398_RAW.tar
curl "${curl_flags[@]}" --continue-at - --output "${ARCHIVE}.partial" "${ARCHIVE_URL}"
test "$(stat -c '%s' "${ARCHIVE}.partial")" -eq 3747420160
mv "${ARCHIVE}.partial" "${ARCHIVE}"
sha256sum "${ARCHIVE}" > "${STAGE}/validation/archive.sha256"

awk -F '\t' '$1 == "File" {print $2}' "${STAGE}/metadata/filelist.txt" \
  | LC_ALL=C sort > "${STAGE}/validation/expected_members.txt"
tar -tf "${ARCHIVE}" | LC_ALL=C sort > "${STAGE}/validation/observed_members.txt"
cmp "${STAGE}/validation/expected_members.txt" "${STAGE}/validation/observed_members.txt"
module load python/3.11.5-GCCcore-13.2.0
python - "${ARCHIVE}" "${STAGE}/metadata/filelist.txt" <<'PY'
from pathlib import Path, PurePosixPath
import csv
import sys
import tarfile

archive, filelist = map(Path, sys.argv[1:])
with filelist.open(encoding="utf-8", newline="") as handle:
    expected = {
        row["Name"]: int(row["Size"])
        for row in csv.DictReader(handle, delimiter="\t")
        if row["#Archive/File"] == "File"
    }
with tarfile.open(archive, mode="r:") as handle:
    members = handle.getmembers()
observed = {}
for member in members:
    path = PurePosixPath(member.name)
    if path.is_absolute() or ".." in path.parts or not member.isfile() or member.name in observed:
        raise SystemExit(f"unsafe or duplicate archive member: {member.name!r}")
    observed[member.name] = member.size
if observed != expected or len(observed) != 30:
    raise SystemExit("archive member census or size differs")
PY

awk -F '\t' '$1 == "File" && $2 ~ /^GSM809/ {print $2}' \
  "${STAGE}/metadata/filelist.txt" > "${STAGE}/validation/human_members.txt"
test "$(wc -l < "${STAGE}/validation/human_members.txt")" -eq 26
tar --extract --file "${ARCHIVE}" --directory "${STAGE}/h5" \
  --no-same-owner --no-same-permissions --verbatim-files-from \
  --files-from "${STAGE}/validation/human_members.txt"
test "$(find "${STAGE}/h5" -maxdepth 1 -type f -name '*.h5' | wc -l)" -eq 26
test "$(find "${STAGE}/h5" -maxdepth 1 -type f -name '*.h5' -printf '%s\n' | awk '{s += $1} END {print s}')" -eq 3392881344
sha256sum "${STAGE}"/h5/*.h5 > "${STAGE}/validation/human_h5.sha256"

module purge
module load h5py/3.9.0-foss-2023a
cd "${ROOT}"
PYTHONPATH=${ROOT}:${PYTHONPATH:-} python -m unittest -v tests.unit.test_audit_gse256398_h5 \
  > "${STAGE}/validation/tests.stdout.txt" 2> "${STAGE}/validation/tests.stderr.txt"
PYTHONPATH=${ROOT}:${PYTHONPATH:-} python "${ROOT}/scripts/audit_gse256398_h5.py" \
  --series-soft "${STAGE}/metadata/GSE256398_family.soft.gz" \
  --filelist "${STAGE}/metadata/filelist.txt" \
  --h5-root "${STAGE}/h5" \
  --output "${STAGE}/audit" \
  > "${STAGE}/validation/audit.stdout.json"
python -m pip freeze --all > "${STAGE}/validation/effective-pip-freeze.txt"
curl --version > "${STAGE}/validation/curl.version.txt"
tar --version > "${STAGE}/validation/tar.version.txt"
sha256sum \
  "${ROOT}/config/datasets/gse256398.toml" \
  "${ROOT}/scripts/audit_gse256398_h5.py" \
  "${ROOT}/tests/unit/test_audit_gse256398_h5.py" \
  "${ROOT}/slurm/acquire_audit_gse256398_h5.sbatch" \
  > "${STAGE}/source.sha256"

module purge
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=${ROOT}/src python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
receipt = json.loads((root / "audit/receipt.json").read_text())
if (
    receipt["status"] != "pass_source_topology_and_h5_structure"
    or receipt["human_donors"] != 26
    or receipt["masld_relevant_donors"] != 17
    or receipt["alcohol_ood_donors"] != 9
    or receipt["human_h5_bytes"] != 3392881344
    or receipt["alcohol_donors_are_masld_negative_controls"]
    or receipt["sealed_or_final_outcomes_read"]
    or receipt["model_training_active"]
):
    raise SystemExit("GSE256398 activation receipt differs")
archive = root / "archive/GSE256398_RAW.tar"
freeze_tree(root, {
    "artifact_class": "gse256398_source_h5_activation",
    "dataset_id": "gse256398",
    "human_donors": 26,
    "masld_relevant_donors": 17,
    "alcohol_ood_donors": 9,
    "biological_unit": "donor",
    "archive_sha256": sha256_file(archive),
    "archive_size_bytes": archive.stat().st_size,
    "native_genome_build": "GRCh38",
    "exact_cell_ranger_reference": "unresolved",
    "external_development_only": True,
    "sealed": False,
    "model_training_active": False,
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE256398 source/H5 activation: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
