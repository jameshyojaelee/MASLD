#!/usr/bin/env bash
#SBATCH --job-name=masld-gse296875-direct
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
FILELIST_URL=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE296nnn/GSE296875/suppl/filelist.txt
EXPECTED_FILELIST_BYTES=2032
EXPECTED_FILELIST_SHA256=983a0764965f832d90ef1e008dcc6bf553450d35c508f4afebf3a1bab2b9768e
EXPECTED_FRAGMENT_BYTES=70290848297
TARGET="${ROOT}/executions/gse296875-fragments-direct-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite direct fragment acquisition: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/fragments" "${STAGE}/validation"
curl --fail --location --silent --show-error \
  --retry 5 --retry-delay 5 --retry-connrefused \
  --output "${STAGE}/filelist.txt" "${FILELIST_URL}"
test "$(stat -c '%s' "${STAGE}/filelist.txt")" -eq "${EXPECTED_FILELIST_BYTES}"
printf '%s  %s\n' "${EXPECTED_FILELIST_SHA256}" "${STAGE}/filelist.txt" \
  | sha256sum --check --strict
awk -F '\t' \
  '$1 == "File" && ($2 ~ /_atac_fragments.tsv.gz$/ || $2 ~ /_atac_fragments.tsv.gz.tbi.gz$/) { print $2 "\t" $4 }' \
  "${STAGE}/filelist.txt" > "${STAGE}/validation/direct_download_roster.tsv"
test "$(wc -l < "${STAGE}/validation/direct_download_roster.tsv")" -eq 16

download_group() {
  local pattern="$1"
  local pids=()
  while IFS=$'\t' read -r name expected_bytes; do
    if [[ ! "${name}" =~ ${pattern} ]]; then
      continue
    fi
    (
      gsm="${name%%_*}"
      url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8979nnn/${gsm}/suppl/${name}"
      partial="${STAGE}/fragments/${name}.partial"
      final="${STAGE}/fragments/${name}"
      complete=false
      for attempt in 1 2 3 4 5; do
        printf 'Downloading %s attempt %s\n' "${name}" "${attempt}"
        if curl --fail --location --silent --show-error \
          --retry 3 --retry-delay 10 --retry-connrefused \
          --continue-at - --output "${partial}" "${url}"; then
          if [[ "$(stat -c '%s' "${partial}")" -eq "${expected_bytes}" ]]; then
            complete=true
            break
          fi
        fi
      done
      if [[ "${complete}" != true ]]; then
        printf 'Direct download incomplete: %s\n' "${name}" >&2
        exit 1
      fi
      mv "${partial}" "${final}"
    ) &
    pids+=("$!")
  done < "${STAGE}/validation/direct_download_roster.tsv"
  for pid in "${pids[@]}"; do
    wait "${pid}"
  done
}

download_group '_atac_fragments\.tsv\.gz$'
download_group '_atac_fragments\.tsv\.gz\.tbi\.gz$'

module load SAMtools/1.23.1
validation_pids=()
for fragment in "${STAGE}"/fragments/*_atac_fragments.tsv.gz; do
  (
    gzip -t "${fragment}"
    bgzip --test "${fragment}"
    compressed_index="${fragment}.tbi.gz"
    gzip -t "${compressed_index}"
    gzip -dc "${compressed_index}" > "${fragment}.tbi"
    tabix -l "${fragment}" > "${fragment}.contigs.txt"
  ) &
  validation_pids+=("$!")
done
for pid in "${validation_pids[@]}"; do
  wait "${pid}"
done

module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}/filelist.txt" "${STAGE}/fragments" \
  "${STAGE}/validation/fragment_headers.json" <<'PY'
from pathlib import Path
import csv
import gzip
import json
import re
import sys

filelist = Path(sys.argv[1]).resolve(strict=True)
root = Path(sys.argv[2]).resolve(strict=True)
output = Path(sys.argv[3])
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
    if path.is_file()
    and not path.is_symlink()
    and (path.name.endswith("_atac_fragments.tsv.gz") or path.name.endswith(".tbi.gz"))
}
if observed != expected:
    raise SystemExit("direct fragment names or sizes differ from GEO filelist")
if sum(size for name, size in observed.items() if name.endswith(".tsv.gz")) != 70_290_848_297:
    raise SystemExit("direct fragment byte total differs")

header_summary = {}
for path in sorted(root.glob("*_atac_fragments.tsv.gz")):
    headers = {}
    primary_contigs = []
    data_rows = 0
    seen_data = False
    previous = None
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                if seen_data:
                    raise SystemExit(f"header after data: {path.name}")
                payload = line[1:].strip()
                if "=" in payload:
                    key, value = payload.split("=", 1)
                    if key == "primary_contig":
                        primary_contigs.append(value)
                    else:
                        headers[key] = value
                continue
            seen_data = True
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise SystemExit(f"fragment data row does not have five fields: {path.name}")
            contig, start_text, end_text, barcode, count_text = fields
            start, end, count = int(start_text), int(end_text), int(count_text)
            if start < 0 or end <= start or count <= 0 or re.fullmatch(r"[ACGT]+-[0-9]+", barcode) is None:
                raise SystemExit(f"invalid fragment data row: {path.name}")
            coordinate = (contig, start, end)
            if previous is not None and contig == previous[0] and coordinate[1:] < previous[1:]:
                raise SystemExit(f"fragment data are not sorted within contig: {path.name}")
            previous = coordinate
            data_rows += 1
            if data_rows == 10_000:
                break
    if data_rows != 10_000:
        raise SystemExit(f"fragment file has fewer than 10,000 data rows: {path.name}")
    if (
        headers.get("pipeline_name") != "cellranger-arc"
        or headers.get("pipeline_version") != "cellranger-arc-2.0.0"
        or headers.get("reference_version") != "2020-A"
        or not {"chr1", "chr22", "chrX", "chrY"}.issubset(primary_contigs)
    ):
        raise SystemExit(f"fragment header contract differs: {path.name}")
    header_summary[path.name] = {
        "pipeline_name": headers["pipeline_name"],
        "pipeline_version": headers["pipeline_version"],
        "reference_version": headers["reference_version"],
        "primary_contigs": primary_contigs,
        "validated_data_rows": data_rows,
    }
output.write_text(json.dumps(header_summary, sort_keys=True, indent=2) + "\n")
PY

sha256sum "${STAGE}"/fragments/* > "${STAGE}/validation/fragment_artifacts.sha256"
cat > "${STAGE}/SOURCE.tsv" <<EOF
field	value
accession	GSE296875
bioproject	PRJNA1261871
source_filelist_url	${FILELIST_URL}
source_filelist_sha256	${EXPECTED_FILELIST_SHA256}
acquisition_mode	direct_GSM_supplement_files
fragment_files	8
fragment_bytes	${EXPECTED_FRAGMENT_BYTES}
source_pipeline	cellranger-arc-2.0.0
source_coordinates	already_adjusted_to_Tn5_cut_sites
fragment_column_5	readSupport_PCR_duplicate_provenance_not_assay_weight
pairing	well_namespaced_same_nucleus_after_frozen_barcode_join
biological_unit	donor
EOF
curl --version > "${STAGE}/validation/curl.version.txt"
bgzip --version > "${STAGE}/validation/bgzip.version.txt"
tabix --version > "${STAGE}/validation/tabix.version.txt"
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
print(
    freeze_tree(
        root,
        {
            "artifact_class": "gse296875_raw_fragment_acquisition",
            "dataset_id": "gse296875",
            "fragment_files": 8,
            "fragment_bytes": 70_290_848_297,
            "acquisition_mode": "direct_GSM_supplement_files",
            "source_pipeline": "cellranger-arc-2.0.0",
            "source_coordinates": "already_adjusted_to_Tn5_cut_sites",
            "automatic_download": False,
            "download_authority": "explicit_user_instruction_to_continue_full_plan_without_further_approval",
            "sealed": False,
        },
    )
)
PY

mv "${STAGE}" "${TARGET}"
printf 'GSE296875 direct fragment acquisition: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
