#!/usr/bin/env bash
#SBATCH --job-name=model-check-953
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
RUNTIME=${ROOT}/executions/environments/transcriptformer-torch2.5.1-21070738
PYTHON=${RUNTIME}/env/bin/python
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

SUBSTRATE=${ROOT}/executions/model-data-950-21181338-showcase-substrate/substrate
SUBSTRATE_SHA=8fed068f5c53d8b0ebbba8c98a734ecd9c99253104a21631e9f58b1c3193bb4b
MANIFEST=${ROOT}/config/showcase_substrate/bayesprism_sources.sha256
RESULTS=${PROJECT}/Analysis/Deconvolution/results
BULK=${PROJECT}/Analysis/Deconvolution/bulk
KDRIVER=${PROJECT}/Analysis/Deconvolution/scripts/21_bayesprism_kallisto.sbatch

TARGET=${ROOT}/executions/model-check-953-${SLURM_JOB_ID}-lineage-variant-overlay
STAGE=${TARGET}.staging
FAILED=${TARGET}.failed

preserve_failed_stage() { if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi; }
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite: %s\n' "${TARGET}" >&2; exit 1
fi
install -d -m 0750 "${STAGE}" "${STAGE}/validation"
cd "${ROOT}"
hostname; "${PYTHON}" -V
bash -n slurm/overlay_lineage_variant_provenance_cpu.sbatch

# 950 and the pinned manifest are READ-ONLY to this step. Checked before, and
# re-checked after, so "additive only" is proven rather than asserted.
cat > "${STAGE}/validation/readonly_before.sha256" <<'SHA'
43137242ad38243387f0b2e85a92f908d62e914a578e7f874669eb8296dbe6ad  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-data-950-21181338-showcase-substrate/ARTIFACTS.json
8fed068f5c53d8b0ebbba8c98a734ecd9c99253104a21631e9f58b1c3193bb4b  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-data-950-21181338-showcase-substrate/substrate/ARTIFACTS.json
1db883f1a5532ccf155bd8fbed54c174904a8429c05ed94b67f4b462f0ccbb95  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/model-data-950-21181338-showcase-substrate/substrate/lineage/receipt.json
8e4ecf0d5729c6f902cead4610077d679892d062a6c22776020470c3e8c7d569  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/config/showcase_substrate/bayesprism_sources.sha256
27c53b301742ef778e40d9b6d1e89a62943c00e13edef5c9479cf14761fd9a11  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Deconvolution/scripts/21_bayesprism_kallisto.sbatch
f22ddd872d681edea9f1c12d5d971c55f66d44d0624d816e317f987ec12c8abb  /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/environments/transcriptformer-torch2.5.1-21070738/ARTIFACTS.json
SHA
sha256sum --check --strict "${STAGE}/validation/readonly_before.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" - <<'PY'
from pathlib import Path
from masld_bench.artifacts import verify_frozen_tree
r = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark")
for t in ("executions/model-data-950-21181338-showcase-substrate",
          "executions/model-data-950-21181338-showcase-substrate/substrate",
          "executions/model-data-950-21181338-showcase-substrate/substrate/lineage",
          "executions/model-data-950-21181338-showcase-substrate/substrate/arm_b"):
    verify_frozen_tree(r / t); print(f"verified {t}")
PY

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" -m unittest -v \
  tests.unit.test_overlay_lineage_variant_provenance \
  tests.unit.test_build_showcase_aspect_lineage_substrate \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

PYTHONPATH="${ROOT}:${ROOT}/src" "${PYTHON}" \
  scripts/overlay_lineage_variant_provenance.py \
  --results-root "${RESULTS}" \
  --bulk-root "${BULK}" \
  --substrate "${SUBSTRATE}" \
  --substrate-artifacts-sha256 "${SUBSTRATE_SHA}" \
  --manifest "${MANIFEST}" \
  --kallisto-driver "${KDRIVER}" \
  --output "${STAGE}/overlay" \
  > "${STAGE}/build_stdout.json" 2> "${STAGE}/build_stderr.txt"

sha256sum --check --strict "${STAGE}/validation/readonly_before.sha256"
cp "${STAGE}/validation/readonly_before.sha256" "${STAGE}/validation/readonly_after.sha256"
printf '950 and the pinned manifest byte-identical before and after: additive only\n'

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt"
sha256sum \
  scripts/overlay_lineage_variant_provenance.py \
  scripts/build_showcase_aspect_lineage_substrate.py \
  tests/unit/test_overlay_lineage_variant_provenance.py \
  slurm/overlay_lineage_variant_provenance_cpu.sbatch \
  "${RUNTIME}/ARTIFACTS.json" \
  > "${STAGE}/source.sha256"

PYTHONPATH="${ROOT}/src" "${PYTHON}" - "${STAGE}" <<'PY'
from pathlib import Path
import json, sys
from masld_bench.artifacts import freeze_tree, verify_frozen_tree

root = Path(sys.argv[1]).resolve(strict=True)
ov = root / "overlay"
r = json.loads((ov / "lineage_variant_overlay.json").read_text())
prov = r["variant_provenance_resolved"]
inv = r["variant_inventory"]
dual = [x for x in inv if x["backup_present"]]

if (
    r["record_is_additive_overlay"] is not True
    or r["bound_substrate"]["substrate_edited"] is not False
    or r["bound_substrate"]["substrate_refrozen"] is not False
    or r["bound_manifest"]["manifest_edited"] is not False
    or len(inv) != 13
    or len(dual) != 5
    or prov["resolved"] is not True
    or prov["signal_1_cohort_membership"]["sets_are_identical"] is not True
    or prov["conclusion"]["star_backup_is_not_a_rounding_variant"] is not True
    or r["crosswalk_finding"]["classification"] != "a missing CHECK, not a missing capability"
    or r["crosswalk_finding"]["no_new_crosswalk_needed"] is not True
    or r["approval"]["user_approved"] is not False
    or any(c["rows_identical"] != 0 for c in r["variant_comparison"])
    or not all(c["roster_identical"] for c in r["variant_comparison"])
):
    raise SystemExit("lineage variant overlay verification differs")

# GSE193066 must be recorded as STAR with no backup: it is the arm a
# both-variant sensitivity cannot cover.
g = next(x for x in inv if x["cohort"] == "GSE193066")
if g["backup_present"] is not False or g.get("canonical_quantification") != "star_human":
    raise SystemExit("GSE193066 variant/quantification record differs")

freeze_tree(ov, {
    "artifact_class": "lineage_variant_provenance_overlay",
    "cohorts": len(inv),
    "cohorts_with_a_backup": len(dual),
    "provenance_resolved": True,
    "substrate_refrozen": False,
    "model_fitted": False,
    "metrics_calculated": False,
    "status": "passed",
})
verify_frozen_tree(ov)
freeze_tree(root, {
    "artifact_class": "lineage_variant_provenance_overlay_campaign",
    "model_fitted": False,
    "metrics_calculated": False,
    "status": "passed",
})
verify_frozen_tree(root)
print(json.dumps({
    "cohorts_with_backup": [x["cohort"] for x in dual],
    "collapse_by_quantification": {
        k: {"cohorts": v["cohorts"], "small_values": v["small_values"]}
        for k, v in r["collapse_structure_is_a_quantification_property"]["by_quantification"].items()},
    "canonical_pairs_with_frozen_axis": r["which_variant_pairs_with_the_frozen_axis"]["frozen_axis_quantification"],
}, indent=2, sort_keys=True))
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'overlay: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json" "${TARGET}/overlay/ARTIFACTS.json"
