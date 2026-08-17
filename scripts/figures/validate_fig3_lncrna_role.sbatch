#!/bin/bash
#SBATCH --job-name=vfig3lnc
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:05:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/vfig3lnc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/vfig3lnc_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FIG3_ROOT="${PROJECT_ROOT}/figures/main/fig3_bulk_transcriptomics"
PDF="${FIG3_ROOT}/panels/fig3_lncrna_role_candidate.pdf"
SOURCE_DIR="${FIG3_ROOT}/source_tables/fig3_lncrna_role"
MANIFEST="${FIG3_ROOT}/manifests/fig3_lncrna_role_candidate_manifest.tsv"

test "$(pdfinfo "${PDF}" | awk '/^Pages:/ {print $2}')" = "1"
strings "${PDF}" | grep -q "/Subtype /TrueType"
! strings "${PDF}" | grep -q "/Subtype /Type3"
test "$(wc -l < "${SOURCE_DIR}/matched_pairs.tsv")" = "369"
awk -F '\t' 'NR > 1 && ($3 < -0.10 || $3 > 0.10) {exit 1}' "${SOURCE_DIR}/match_balance.tsv"

python3 - "${MANIFEST}" <<'PY'
import csv
import hashlib
import pathlib
import sys

with pathlib.Path(sys.argv[1]).open(newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        path = pathlib.Path(row["path"])
        assert path.stat().st_size == int(row["size_bytes"]), path
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"], path
PY

python3 scripts/manuscript/validate_resource_scope.py
echo "FIG3_LNCRNA_ROLE_VALIDATION PASS"
