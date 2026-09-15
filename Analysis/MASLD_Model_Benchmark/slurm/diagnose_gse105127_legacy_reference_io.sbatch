#!/usr/bin/env bash
#SBATCH --job-name=data-check-211
#SBATCH --partition=io
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/gse105127-production-9675cd48993468d4/attempts/reference-21083577.failed/output.staging/human_g1k_v37.fasta.gz
SOURCE_SHA256=8b6c538abf0dd92d3f3020f36cc1dd67ce004ffa421c2781205f1eb690bdb442
TARGET=${ROOT}/executions/data-check-211-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() {
  status=$?
  if [[ ${status} -ne 0 && -d ${STAGE} && ! -e ${FAILED} ]]; then
    mv "${STAGE}" "${FAILED}"
  fi
  exit "${status}"
}
trap preserve_failed_stage EXIT
if [[ -e ${TARGET} || -e ${STAGE} || -e ${FAILED} ]]; then
  printf 'Refusing to overwrite legacy reference diagnostic.\n' >&2
  exit 1
fi
install -d -m 0750 "${STAGE}"
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}" | sha256sum --check --strict
if [[ $(stat -c '%s' "${SOURCE}") -ne 892331003 ]]; then
  printf 'Legacy compressed reference size differs.\n' >&2
  exit 1
fi

set +e
gzip -cd "${SOURCE}" 2>"${STAGE}/gzip.stderr.txt" \
  | python3 -c 'import hashlib,json,sys
h1=hashlib.md5(); h2=hashlib.sha256(); total=0
for block in iter(lambda: sys.stdin.buffer.read(8*1024*1024), b""):
    total += len(block); h1.update(block); h2.update(block)
with open(sys.argv[1], "x", encoding="utf-8") as handle:
    json.dump({"bytes": total, "md5": h1.hexdigest(), "sha256": h2.hexdigest()}, handle, sort_keys=True)
    handle.write("\n")' "${STAGE}/uncompressed_stream.json"
statuses=("${PIPESTATUS[@]}")
set -e
if [[ ${statuses[0]} -ne 2 || ${statuses[1]} -ne 0 ]]; then
  printf 'Legacy gzip/Python stream statuses differ: %s %s\n' "${statuses[0]}" "${statuses[1]}" >&2
  exit 1
fi
if ! grep -Fq 'trailing garbage ignored' "${STAGE}/gzip.stderr.txt"; then
  printf 'Expected legacy trailing-payload warning is absent.\n' >&2
  exit 1
fi
if grep -Eqi 'unexpected end|crc error|invalid compressed data|format violated' "${STAGE}/gzip.stderr.txt"; then
  printf 'Legacy stream has a fatal gzip diagnostic.\n' >&2
  exit 1
fi

module purge
module load python/3.11.5-GCCcore-13.2.0
cd "${ROOT}"
PYTHONPATH=src python - "${STAGE}" "${SOURCE_SHA256}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree, write_json_exclusive

stage = Path(sys.argv[1]).resolve(strict=True)
stream = json.loads((stage / "uncompressed_stream.json").read_text(encoding="utf-8"))
if stream.get("md5") != "0ce84c872fc0072a885926823dcd0338" or stream.get("bytes", 0) < 3_000_000_000:
    raise SystemExit("canonical uncompressed FASTA identity differs")
write_json_exclusive(stage / "diagnostic_receipt.json", {
    "schema_version": "masld-bench-gse105127-legacy-reference-stream-diagnostic-v1",
    "status": "pass_exact_legacy_stream",
    "compressed_bytes": 892_331_003,
    "compressed_sha256": sys.argv[2],
    "gzip_exit_status": 2,
    "legacy_warning": "trailing garbage ignored",
    "uncompressed_bytes": stream["bytes"],
    "uncompressed_md5": stream["md5"],
    "uncompressed_sha256": stream["sha256"],
    "fatal_gzip_diagnostic": False,
    "fit_or_score_performed": False,
    "labels_accessed": False,
})
freeze_tree(stage, {
    "artifact_class": "gse105127_legacy_reference_stream_diagnostic",
    "compressed_sha256": sys.argv[2],
    "uncompressed_md5": stream["md5"],
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
