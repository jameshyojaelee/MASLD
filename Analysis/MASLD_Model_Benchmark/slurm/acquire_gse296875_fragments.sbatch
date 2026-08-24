#!/usr/bin/env bash
#SBATCH --job-name=masld-gse296875-fragments
#SBATCH --partition=io
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
URL=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296875/suppl/GSE296875_RAW.tar
FILELIST_URL=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296875/suppl/filelist.txt
EXPECTED_ARCHIVE_BYTES=72391659520
EXPECTED_FILELIST_BYTES=2032
EXPECTED_FILELIST_SHA256=983a0764965f832d90ef1e008dcc6bf553450d35c508f4afebf3a1bab2b9768e
TARGET="${ROOT}/executions/gse296875-fragments-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
ARCHIVE="${STAGE}/archive/GSE296875_RAW.tar"
PARTIAL="${ARCHIVE}.partial"

cleanup() {
  if [[ -d "${STAGE}" ]]; then
    rm -rf -- "${STAGE}"
  fi
}
trap cleanup EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite acquisition output: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/archive" "${STAGE}/fragments" "${STAGE}/validation"

curl --fail --location --silent --show-error \
  --retry 5 --retry-delay 5 --retry-connrefused \
  "${FILELIST_URL}" \
  --output "${STAGE}/filelist.txt"
test "$(stat -c '%s' "${STAGE}/filelist.txt")" -eq "${EXPECTED_FILELIST_BYTES}"
printf '%s  %s\n' "${EXPECTED_FILELIST_SHA256}" "${STAGE}/filelist.txt" \
  | sha256sum --check --strict

download_complete=false
for attempt in 1 2 3 4 5; do
  printf 'Archive download attempt %s\n' "${attempt}"
  if curl --fail --location --show-error \
    --retry 3 --retry-delay 10 --retry-connrefused \
    --continue-at - \
    "${URL}" \
    --output "${PARTIAL}"; then
    if [[ "$(stat -c '%s' "${PARTIAL}")" -eq "${EXPECTED_ARCHIVE_BYTES}" ]]; then
      download_complete=true
      break
    fi
  fi
done
if [[ "${download_complete}" != true ]]; then
  printf 'Archive download did not reach expected size %s; observed %s\n' \
    "${EXPECTED_ARCHIVE_BYTES}" "$(stat -c '%s' "${PARTIAL}" 2>/dev/null || printf 0)" >&2
  exit 1
fi
mv "${PARTIAL}" "${ARCHIVE}"
sha256sum "${ARCHIVE}" > "${STAGE}/validation/archive.sha256"

awk -F '\t' '$1 == "File" { print $2 }' "${STAGE}/filelist.txt" \
  | LC_ALL=C sort > "${STAGE}/validation/expected_archive_members.txt"
tar -tf "${ARCHIVE}" | LC_ALL=C sort \
  > "${STAGE}/validation/observed_archive_members.txt"
cmp \
  "${STAGE}/validation/expected_archive_members.txt" \
  "${STAGE}/validation/observed_archive_members.txt"

module load python/3.11.5-GCCcore-13.2.0
python - "${ARCHIVE}" "${STAGE}/filelist.txt" <<'PY'
from pathlib import Path, PurePosixPath
import csv
import sys
import tarfile

archive = Path(sys.argv[1]).resolve(strict=True)
filelist = Path(sys.argv[2]).resolve(strict=True)
expected = {}
with filelist.open(encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["#Archive/File"] == "File":
            expected[row["Name"]] = int(row["Size"])
with tarfile.open(archive, mode="r:") as handle:
    members = handle.getmembers()
observed = {}
for member in members:
    path = PurePosixPath(member.name)
    if path.is_absolute() or ".." in path.parts or not member.isfile():
        raise SystemExit(f"unsafe archive member: {member.name!r}")
    if member.name in observed:
        raise SystemExit(f"duplicate archive member: {member.name!r}")
    observed[member.name] = member.size
if observed != expected:
    raise SystemExit("archive member sizes differ from authoritative filelist")
if len(observed) != 24:
    raise SystemExit(f"expected 24 regular archive members; observed {len(observed)}")
PY

awk -F '\t' \
  '$1 == "File" && ($2 ~ /_atac_fragments.tsv.gz$/ || $2 ~ /_atac_fragments.tsv.gz.tbi.gz$/) { print $2 }' \
  "${STAGE}/filelist.txt" \
  | LC_ALL=C sort > "${STAGE}/validation/fragment_archive_members.txt"
test "$(wc -l < "${STAGE}/validation/fragment_archive_members.txt")" -eq 16
tar \
  --extract \
  --file "${ARCHIVE}" \
  --directory "${STAGE}/fragments" \
  --no-same-owner \
  --no-same-permissions \
  --verbatim-files-from \
  --files-from "${STAGE}/validation/fragment_archive_members.txt"

python - "${STAGE}/filelist.txt" "${STAGE}/fragments" <<'PY'
from pathlib import Path
import csv
import sys

filelist = Path(sys.argv[1]).resolve(strict=True)
root = Path(sys.argv[2]).resolve(strict=True)
expected = {}
with filelist.open(encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        name = row["Name"]
        if name.endswith("_atac_fragments.tsv.gz") or name.endswith(
            "_atac_fragments.tsv.gz.tbi.gz"
        ):
            expected[name] = int(row["Size"])
observed = {
    path.name: path.stat().st_size
    for path in root.iterdir()
    if path.is_file() and not path.is_symlink()
}
if observed != expected:
    raise SystemExit("extracted fragment file names or sizes differ")
if sum(size for name, size in observed.items() if name.endswith(".tsv.gz")) != 70290848297:
    raise SystemExit("eight-fragment byte total differs")
PY

module load SAMtools/1.23.1
for fragment in "${STAGE}"/fragments/*_atac_fragments.tsv.gz; do
  gzip -t "${fragment}"
  bgzip --test "${fragment}"
  compressed_index="${fragment}.tbi.gz"
  gzip -t "${compressed_index}"
  gzip -dc "${compressed_index}" > "${fragment}.tbi"
  tabix -l "${fragment}" > "${fragment}.contigs.txt"
  for required in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 \
    chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX; do
    if ! awk -v contig="${required}" '$0 == contig { found = 1 } END { exit !found }' \
      "${fragment}.contigs.txt"; then
      printf 'Required contig missing from %s index: %s\n' "${fragment}" "${required}" >&2
      exit 1
    fi
  done
done

python - "${STAGE}/fragments" <<'PY'
from pathlib import Path
import gzip
import sys

root = Path(sys.argv[1]).resolve(strict=True)
for path in sorted(root.glob("*_atac_fragments.tsv.gz")):
    checked = 0
    previous = None
    with gzip.open(path, mode="rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise SystemExit(f"fragment row has fewer than five fields: {path.name}")
            chromosome, start_text, end_text, barcode, count_text = fields[:5]
            start, end, count = int(start_text), int(end_text), int(count_text)
            if not chromosome.startswith("chr") or start < 0 or end <= start or count <= 0:
                raise SystemExit(f"invalid fragment row: {path.name}")
            key = (chromosome, start, end)
            if previous is not None and key < previous:
                raise SystemExit(f"fragment sample is not coordinate sorted: {path.name}")
            if not barcode or barcode.startswith("well"):
                raise SystemExit(f"fragment barcode is not raw well-local identity: {path.name}")
            previous = key
            checked += 1
            if checked == 10000:
                break
    if checked != 10000:
        raise SystemExit(f"fragment file has fewer than 10,000 rows: {path.name}")
PY

sha256sum "${STAGE}"/fragments/* \
  > "${STAGE}/validation/fragment_artifacts.sha256"
cat > "${STAGE}/SOURCE.tsv" <<EOF
field	value
accession	GSE296875
bioproject	PRJNA1261871
source_url	${URL}
source_archive_expected_bytes	${EXPECTED_ARCHIVE_BYTES}
source_filelist_url	${FILELIST_URL}
source_filelist_sha256	${EXPECTED_FILELIST_SHA256}
fragment_files	8
fragment_bytes	70290848297
pairing	well_namespaced_same_nucleus_after_frozen_barcode_join
biological_unit	donor
EOF
curl --version > "${STAGE}/validation/curl.version.txt"
tar --version > "${STAGE}/validation/tar.version.txt"
bgzip --version > "${STAGE}/validation/bgzip.version.txt"
tabix --version > "${STAGE}/validation/tabix.version.txt"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree, sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
archive = root / "archive" / "GSE296875_RAW.tar"
freeze_tree(
    root,
    {
        "artifact_class": "gse296875_raw_fragment_acquisition",
        "dataset_id": "gse296875",
        "archive_sha256": sha256_file(archive),
        "archive_size_bytes": archive.stat().st_size,
        "fragment_files": 8,
        "fragment_bytes": 70290848297,
        "automatic_download": False,
        "download_authority": "explicit_user_instruction_to_continue_full_plan_without_further_approval",
        "sealed": False,
    },
)
PY

mv "${STAGE}" "${TARGET}"
trap - EXIT
printf 'GSE296875 fragment acquisition: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
