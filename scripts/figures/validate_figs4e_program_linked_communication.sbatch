#!/bin/bash
#SBATCH --job-name=validate_figs4e_ccc
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/validate_figs4e_ccc_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/validate_figs4e_ccc_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FIGURE_ROOT="${PROJECT_ROOT}/figures/main/fig4_singlecell_programs"
MANIFEST="${FIGURE_ROOT}/manifests/current_candidate_manifest.tsv"
SUPP_PDF="${FIGURE_ROOT}/panels/supplementary/figs4e_program_linked_communication.pdf"
SOURCE_CANDIDATE="${PROJECT_ROOT}/figures/candidates/fig4cd-ccc-bridge-2026-08-17-v13"

test ! -e "${FIGURE_ROOT}/panels/fig4d_communication.pdf"
test -s "${FIGURE_ROOT}/panels/_superseded_main_communication_2026-08-17_v1/fig4d_communication.pdf"
test -s "${SUPP_PDF}"
test "$(sha256sum "${SUPP_PDF}" | awk '{print $1}')" = \
  "$(sha256sum "${SOURCE_CANDIDATE}/panels/fig4d_program_linked_communication.pdf" | awk '{print $1}')"

while IFS=$'\t' read -r role relative_path expected_size expected_sha state; do
  if [[ "${role}" == "role" ]]; then
    continue
  fi
  artifact="${FIGURE_ROOT}/${relative_path}"
  test -s "${artifact}"
  test "$(stat -c %s "${artifact}")" = "${expected_size}"
  test "$(sha256sum "${artifact}" | awk '{print $1}')" = "${expected_sha}"
done < "${MANIFEST}"

test "$(pdfinfo "${SUPP_PDF}" | awk '/^Pages:/ {print $2}')" = "1"
gs -q -dNOPAUSE -dBATCH -dPDFDEBUG -sDEVICE=nullpage "${SUPP_PDF}"

python3 - "${SUPP_PDF}" "${FIGURE_ROOT}" <<'PY'
import csv
from pathlib import Path
import re
import sys

pdf = Path(sys.argv[1])
root = Path(sys.argv[2])
raw = pdf.read_bytes()
assert re.search(rb"/Subtype\s*/Type3\b", raw) is None
assert re.search(rb"/FontType\s+3\b", raw) is None

source = root / "source_tables/current_candidate/figs4e_program_linked_communication"
with (source / "fig4d_program_ccc_associations.tsv").open(newline="") as handle:
    associations = list(csv.DictReader(handle, delimiter="\t"))
assert len(associations) == 26
assert sum(row["testable"] == "TRUE" for row in associations) == 24
assert sum(row["evidence_state"] == "supported" for row in associations) == 1

with (source / "fig4d_supported_bridge_partial_residuals.tsv").open(newline="") as handle:
    partial = list(csv.DictReader(handle, delimiter="\t"))
assert len(partial) == 60
assert len({row["donor"] for row in partial}) == 60
PY

python3 scripts/manuscript/validate_resource_scope.py
echo "FIGURE_S4E_ROUTING_VALIDATION PASS"
