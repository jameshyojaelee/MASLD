#!/usr/bin/env bash
#SBATCH --job-name=gse267145-coord-freeze
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
SUPPLEMENT="${ROOT}/executions/gse267145-supplement3-21070402"
METADATA="${ROOT}/executions/geo-multimodal-source-audit-21064310"
MATRIX="${ROOT}/executions/gse267145-matrix-audit-21064333"
JOIN="${ROOT}/executions/gse267145-authoritative-join-21064930"
PRIOR="${ROOT}/executions/gse267145-coordinate-activation-21070593"
CHROMHMM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/roadmap_liver_E066/E066_15_chromHMM_hg38.bed.gz
TARGET="${ROOT}/executions/gse267145-coordinate-convention-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite coordinate-convention freeze: %s\n' "${TARGET}" >&2; exit 1
fi
cd "${ROOT}"
module load python/3.11.5-GCCcore-13.2.0
mkdir -p "${STAGE}/validation"
# The full module carries a PRE-EXISTING failure unrelated to this freeze:
# test_declared_task_is_development_only_and_fail_closed raises "promotion gate SHA-256
# differs" because config/evaluation/paired_bulk_rna_h3k27ac_promotion_gate.json no longer
# hashes to the value pinned in its TaskSpec. That is config drift, recorded for the owner,
# not fixed here. The freeze runs the coordinate-claims tests and the structural-evidence
# tests, and records the full-module result beside them without gating on it.
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse267145_coordinate_activation.GSE267145CoordinateActivationTests.test_bed_mention_and_in_bounds_do_not_resolve_semantics \
  tests.unit.test_audit_gse267145_coordinate_activation.GSE267145CoordinateActivationTests.test_only_direct_coordinate_claim_is_recognized \
  tests.unit.test_audit_gse267145_coordinate_activation.GSE267145CoordinateStructuralEvidenceTests \
  > "${STAGE}/validation/unit_tests.txt" 2>&1
PYTHONPATH="src:." python -m unittest -v \
  tests.unit.test_audit_gse267145_coordinate_activation \
  > "${STAGE}/validation/unit_tests_full_module_informational.txt" 2>&1 || \
  printf 'full module has a pre-existing failure; see validation/unit_tests_full_module_informational.txt\n' >&2
PYTHONPATH="src" python - "${SUPPLEMENT}" "${METADATA}" "${MATRIX}" "${JOIN}" "${PRIOR}" <<'PY'
from pathlib import Path
import sys
from masld_bench.artifacts import verify_frozen_tree
for value in sys.argv[1:]:
    verify_frozen_tree(Path(value).resolve(strict=True))
print("input frozen-tree verification passed")
PY
sha256sum "${CHROMHMM}" > "${STAGE}/validation/chromhmm_input.sha256"
PYTHONPATH="src:." python scripts/freeze_gse267145_coordinate_convention.py \
  --workbook "${SUPPLEMENT}/raw/media-3.xlsx" \
  --matrix "${MATRIX}/raw/GSE267119_H3K27ac.txt.gz" \
  --chromhmm-bed "${CHROMHMM}" \
  --article-xml "${JOIN}/raw/PMC11071482.xml" \
  --geo-soft "${METADATA}/raw/GSE267119_family.soft.gz" \
  --prior-audit "${PRIOR}/coordinate_audit.json" \
  --output "${STAGE}/freeze" | tee "${STAGE}/validation/freeze_stdout.json"
python -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true
sha256sum \
  scripts/audit_gse267145_coordinate_activation.py \
  scripts/freeze_gse267145_coordinate_convention.py \
  tests/unit/test_audit_gse267145_coordinate_activation.py \
  slurm/freeze_gse267145_coordinate_convention.sbatch \
  > "${STAGE}/source.sha256"
PYTHONPATH="src" python - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree
root = Path(sys.argv[1]).resolve(strict=True)
fz = json.loads((root / "freeze/convention_freeze.json").read_text(encoding="utf-8"))
if not fz["frozen"]:
    raise SystemExit("structural tests did not agree; convention stays unresolved; nothing frozen")
for flag in ("coordinate_semantics_resolved_by_primary_source",
             "bounds_compatibility_used_as_semantic_evidence",
             "tool_default_or_file_extension_used_as_semantic_evidence"):
    if fz[flag] is not False:
        raise SystemExit(f"{flag} must be False")
if fz["coordinate_semantics_resolved_flag_in_prior_audit"] is not False:
    raise SystemExit("prior audit flag must remain False")
freeze_tree(root, {
    "artifact_class": "gse267145_coordinate_convention_structural_freeze",
    "cohort_family_id": "gse267145_znf469_human_liver",
    "coordinate_convention_frozen": fz["coordinate_convention_frozen"],
    "resolution_class": fz["resolution_class"],
    "primary_source_statement_found": False,
    "coordinate_semantics_resolved": False,
    "tier1_gene_and_interval_operations_allowed": True,
    "tier2_called_peak_base_resolution_allowed": False,
    "champion_eligible": False,
    "status": "passed",
})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE267145 coordinate convention freeze: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
