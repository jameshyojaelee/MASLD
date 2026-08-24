#!/usr/bin/env bash
#SBATCH --job-name=masld-encode-ccres
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
ACCESSION=ENCFF420VPZ
URL=https://www.encodeproject.org/files/ENCFF420VPZ/@@download/ENCFF420VPZ.bed.gz
METADATA_URL=https://www.encodeproject.org/files/ENCFF420VPZ/?format=json
EXPECTED_BYTES=32795453
EXPECTED_MD5=ea1820bd71581ba01c94fc600fbf794d
TARGET="${ROOT}/executions/encode-v4-ccres-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
BED="${STAGE}/${ACCESSION}.bed.gz"

if [[ -e "${TARGET}" || -e "${STAGE}" ]]; then
  printf 'Refusing to overwrite ENCODE cCRE acquisition: %s\n' "${TARGET}" >&2
  exit 1
fi
install -d -m 0750 "${STAGE}/validation"
curl --fail --location --silent --show-error \
  --retry 5 --retry-delay 5 --retry-connrefused \
  --output "${BED}" "${URL}"
curl --fail --location --silent --show-error \
  --retry 5 --retry-delay 5 --retry-connrefused \
  --header 'Accept: application/json' \
  --output "${STAGE}/${ACCESSION}.metadata.json" "${METADATA_URL}"
test "$(stat -c '%s' "${BED}")" -eq "${EXPECTED_BYTES}"
printf '%s  %s\n' "${EXPECTED_MD5}" "${BED}" | md5sum --check --strict
gzip -t "${BED}"

module load python/3.11.5-GCCcore-13.2.0
python - "${BED}" "${STAGE}/${ACCESSION}.metadata.json" \
  "${STAGE}/validation/content.json" <<'PY'
from collections import Counter
import gzip
import json
from pathlib import Path
import sys

bed = Path(sys.argv[1]).resolve(strict=True)
metadata_path = Path(sys.argv[2]).resolve(strict=True)
output = Path(sys.argv[3])
metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
required = {
    "accession": "ENCFF420VPZ",
    "assembly": "GRCh38",
    "file_format": "bed",
    "file_size": 32_795_453,
    "md5sum": "ea1820bd71581ba01c94fc600fbf794d",
    "output_type": "candidate Cis-Regulatory Elements",
    "status": "released",
}
for key, value in required.items():
    if metadata.get(key) != value:
        raise SystemExit(f"ENCODE metadata differs for {key}: {metadata.get(key)!r}")
if metadata.get("superseded_by") not in ([], None):
    raise SystemExit("ENCODE cCRE file has been superseded")
if "ENCODE v4" not in metadata.get("encyclopedia_version", []):
    raise SystemExit("ENCODE encyclopedia version differs")

contigs = Counter()
classes = Counter()
identifiers = set()
seen_contigs = set()
current_contig = None
previous_coordinate = None
rows = 0
with gzip.open(bed, "rt", encoding="utf-8") as handle:
    for line_number, line in enumerate(handle, start=1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            raise SystemExit(f"cCRE row {line_number} has fewer than BED9 fields")
        contig, start_text, end_text, identifier = fields[:4]
        start, end = int(start_text), int(end_text)
        if start < 0 or end <= start or not identifier:
            raise SystemExit(f"cCRE row {line_number} has invalid coordinates or ID")
        if identifier in identifiers:
            raise SystemExit(f"duplicate cCRE ID: {identifier}")
        identifiers.add(identifier)
        if contig != current_contig:
            if contig in seen_contigs:
                raise SystemExit("ENCODE cCRE contig appears in multiple blocks")
            seen_contigs.add(contig)
            current_contig = contig
            previous_coordinate = None
        coordinate = (start, end)
        if previous_coordinate is not None and coordinate < previous_coordinate:
            raise SystemExit("ENCODE cCRE BED is not sorted within contig")
        previous_coordinate = coordinate
        contigs[contig] += 1
        classes[fields[8]] += 1
        rows += 1
if rows < 1_000_000:
    raise SystemExit(f"ENCODE cCRE row count is unexpectedly small: {rows}")
output.write_text(
    json.dumps(
        {
            "rows": rows,
            "unique_ids": len(identifiers),
            "contigs": dict(sorted(contigs.items())),
            "bed9_field_9_classes": dict(sorted(classes.items())),
        },
        sort_keys=True,
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)
PY

sha256sum "${BED}" "${STAGE}/${ACCESSION}.metadata.json" \
  > "${STAGE}/validation/source.sha256"
curl --version > "${STAGE}/validation/curl.version.txt"
cat > "${STAGE}/SOURCE.tsv" <<EOF
field	value
accession	${ACCESSION}
assembly	GRCh38
encyclopedia_version	ENCODE_v4
source_url	${URL}
metadata_url	${METADATA_URL}
expected_bytes	${EXPECTED_BYTES}
expected_md5	${EXPECTED_MD5}
status	released
rights	unrestricted_use_per_ENCODE_data_use_policy
rights_url	https://www.encodeproject.org/about/data-use-policy/
benchmark_role	outcome_independent_held_chromosome_evaluation_reference
EOF
PYTHONPATH="${ROOT}/src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys

from masld_bench.artifacts import freeze_tree
from masld_bench.hashing import sha256_file

root = Path(sys.argv[1]).resolve(strict=True)
content = json.loads((root / "validation" / "content.json").read_text())
print(
    freeze_tree(
        root,
        {
            "artifact_class": "outcome_independent_genomic_reference",
            "accession": "ENCFF420VPZ",
            "assembly": "GRCh38",
            "encyclopedia_version": "ENCODE_v4",
            "rows": content["rows"],
            "source_sha256": sha256_file(root / "ENCFF420VPZ.bed.gz"),
            "sealed_outcome_used": False,
        },
    )
)
PY

mv "${STAGE}" "${TARGET}"
printf 'ENCODE v4 cCRE acquisition: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
