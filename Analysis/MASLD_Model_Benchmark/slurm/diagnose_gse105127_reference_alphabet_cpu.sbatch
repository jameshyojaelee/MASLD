#!/usr/bin/env bash
#SBATCH --job-name=data-check-216
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=executions/logs/data-check-216-%j.out
#SBATCH --error=executions/logs/data-check-216-%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SOURCE=${ROOT}/executions/gse105127-production-9675cd48993468d4/attempts/reference-exact-21086153.failed/output.staging/human_g1k_v37.fasta
SOURCE_SHA256=2f9cd9e853a9284c53884e6a551b1c7284795dd053f255d630aeeb114d1fa81f
TARGET=${ROOT}/executions/data-check-216-${SLURM_JOB_ID}
STAGE=${TARGET}.staging
cd "${ROOT}"
if [[ -e ${TARGET} || -e ${STAGE} ]]; then
  printf 'Refusing to overwrite FASTA alphabet diagnostic.\n' >&2
  exit 1
fi
mkdir -p "${STAGE}"
printf '%s  %s\n' "${SOURCE_SHA256}" "${SOURCE}" | sha256sum --check --strict
PYTHONPATH="${ROOT}/src" \
  "${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python" \
  - "${SOURCE}" "${STAGE}" <<'PY'
from collections import Counter
import json
from pathlib import Path
import sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive

source = Path(sys.argv[1]).resolve(strict=True)
stage = Path(sys.argv[2]).resolve(strict=True)
counts = Counter()
sequence_rows = 0
header_rows = 0
with source.open("rb") as handle:
    for line in handle:
        if line.startswith(b">"):
            header_rows += 1
            continue
        sequence = line.strip()
        sequence_rows += 1
        counts.update(sequence)
observed = {
    chr(value): count for value, count in sorted(counts.items())
}
write_json_exclusive(stage / "alphabet_receipt.json", {
    "schema_version": "masld-bench-gse105127-reference-alphabet-diagnostic-v1",
    "status": "passed_outcome_free_census",
    "source_bytes": source.stat().st_size,
    "source_sha256": "2f9cd9e853a9284c53884e6a551b1c7284795dd053f255d630aeeb114d1fa81f",
    "header_rows": header_rows,
    "sequence_rows": sequence_rows,
    "alphabet_counts": observed,
    "alphabet": sorted(observed),
    "labels_accessed": False,
    "fit_or_score_performed": False,
})
freeze_tree(stage, {
    "artifact_class": "gse105127_reference_alphabet_diagnostic",
    "labels_accessed": False,
    "fit_or_score_performed": False,
    "status": "passed",
})
verify_frozen_tree(stage)
PY
mv "${STAGE}" "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
