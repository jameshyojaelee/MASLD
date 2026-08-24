#!/usr/bin/env bash
#SBATCH --job-name=masld-encode-blacklist
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
ACCESSION=ENCFF356LFX
TARGET="${ROOT}/executions/encode-grch38-blacklist-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
METADATA_URL="https://www.encodeproject.org/files/${ACCESSION}/?format=json"
DOWNLOAD_URL="https://www.encodeproject.org/files/${ACCESSION}/@@download/${ACCESSION}.bed.gz"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite ENCODE blacklist acquisition: %s\n' "${TARGET}" >&2
  exit 1
fi
mkdir -m 750 "${STAGE}"
module load curl/8.7.1-GCCcore-13.3.0
curl --fail --location --retry 4 --retry-delay 5 \
  --header 'Accept: application/json' \
  "${METADATA_URL}" --output "${STAGE}/${ACCESSION}.metadata.json"
curl --fail --location --retry 4 --retry-delay 5 \
  "${DOWNLOAD_URL}" --output "${STAGE}/${ACCESSION}.bed.gz"
printf '%s  %s\n' 393688b4f06c9ce26165d47433dd8c37 \
  "${STAGE}/${ACCESSION}.bed.gz" | md5sum --check --strict
gzip -t "${STAGE}/${ACCESSION}.bed.gz"

module load python/3.11.5-GCCcore-13.2.0
python - "${STAGE}" "${ACCESSION}" <<'PY'
from pathlib import Path
import gzip
import json
import sys

root = Path(sys.argv[1]).resolve(strict=True)
accession = sys.argv[2]
metadata = json.loads((root / f"{accession}.metadata.json").read_text())
if (
    metadata.get("accession") != accession
    or metadata.get("assembly") != "GRCh38"
    or metadata.get("status") != "released"
    or metadata.get("file_format") != "bed"
    or metadata.get("file_format_type") != "bed3"
    or metadata.get("output_type") != "exclusion list regions"
    or metadata.get("md5sum") != "393688b4f06c9ce26165d47433dd8c37"
):
    raise SystemExit("ENCODE blacklist metadata differs from the contract")
rows = 0
previous = None
with gzip.open(root / f"{accession}.bed.gz", "rt", encoding="utf-8") as handle:
    for line_number, line in enumerate(handle, start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 3:
            raise SystemExit(f"blacklist row {line_number} is not BED3")
        contig, start_text, end_text = fields
        start, end = int(start_text), int(end_text)
        if start < 0 or end <= start:
            raise SystemExit(f"blacklist row {line_number} is invalid")
        key = (contig, start, end)
        if previous is not None and key == previous:
            raise SystemExit("blacklist contains an exact duplicate row")
        previous = key
        rows += 1
if rows == 0:
    raise SystemExit("blacklist is empty")
(root / "validation.json").write_text(
    json.dumps({"accession": accession, "assembly": "GRCh38", "rows": rows}, sort_keys=True, indent=2) + "\n"
)
PY
curl --version > "${STAGE}/curl.version.txt"
python -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum "${ROOT}/slurm/acquire_encode_grch38_blacklist.sbatch" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
validation = json.loads((root / "validation.json").read_text())
print(
    freeze_tree(
        root,
        {
            "artifact_class": "encode_grch38_blacklist",
            "accession": validation["accession"],
            "assembly": validation["assembly"],
            "rows": validation["rows"],
            "license": "ENCODE_unrestricted_data_use_policy",
        },
    )
)
PY
mv "${STAGE}" "${TARGET}"
printf 'ENCODE GRCh38 blacklist: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
