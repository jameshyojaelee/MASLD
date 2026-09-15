#!/usr/bin/env bash
#SBATCH --job-name=model-data-002
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
LOCK="${PROJECT}/Analysis/SingleCell/continual_integration/reference/gse296875_source_lock_v44.json"
CORGI="${ROOT}/executions/corgi-gse296875-context-counts-21066278"
SPEC="${ROOT}/executions/gse296875-phenotype-campaign-spec-21107179"
TARGET="${ROOT}/executions/gse296875-donor-lineage-pseudobulk-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite pseudobulk artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
mkdir -p "${STAGE}/validation"

# Fail closed on every frozen input, including the spec frozen by model-check-001.
cat > "${STAGE}/validation/bound_fixtures.sha256" <<'SUMS'
582387b78ec25215d45820c24aace12f130b826e91755bf19ac4581730611e64  executions/gse296875-fragment-membership-21063829/ARTIFACTS.json
ead9a54099fe2dfb47517497c072feec53fa5ebf9c0392e38ee34a164f8d64e8  executions/gse296875-barcode-donor-join-20260825/ARTIFACTS.json
6518494a7249cb1eac81fbad75989e7dfe0269b611254d5097736b9f5b1da140  executions/gse296875-phenotype-endpoints-20260825/ARTIFACTS.json
e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416  executions/gse296875-donor-folds-20260825/ARTIFACTS.json
06986191ab76f9d8bfba0dbda109ecfdeff05354f37a34fb9a5d4bf49dd50fae  executions/gse296875-donor-supplement-20260825-v2/ARTIFACTS.json
6946c329e74cd1a6cd72fe37faaf805183ff24c729ef9e549da11c49ca02fef8  executions/corgi-gse296875-context-counts-21066278/ARTIFACTS.json
474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58  executions/gse296875-phenotype-campaign-spec-21107179/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

"${PYTHON}" -m unittest -v tests.unit.test_build_gse296875_donor_lineage_pseudobulk \
  > "${STAGE}/validation/unit_tests.txt" 2>&1

"${PYTHON}" scripts/build_gse296875_donor_lineage_pseudobulk.py \
  --root "${ROOT}" \
  --source-lock "${LOCK}" \
  --corgi-root "${CORGI}" \
  --output "${STAGE}/pseudobulk" \
  > "${STAGE}/validation/build_stdout.json"

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true

sha256sum \
  scripts/build_gse296875_donor_lineage_pseudobulk.py \
  tests/unit/test_build_gse296875_donor_lineage_pseudobulk.py \
  slurm/build_gse296875_donor_lineage_pseudobulk_cpu.sbatch \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python3 - "${STAGE}" "${SPEC}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
spec = json.loads(
    (Path(sys.argv[2]) / "spec" / "campaign_spec.json").read_text()
)
audit = json.loads((root / "pseudobulk" / "audit.json").read_text())

if audit["outcomes_read"] or audit["model_training_activated"]:
    raise SystemExit("pseudobulk receipt claims an outcome was read")
if audit["primary_lineages"] != spec["lineages"]["primary"]:
    raise SystemExit("pseudobulk lineage partition differs from the frozen spec")
if audit["minimum_nuclei_per_unit"] != spec["fixture"]["minimum_nuclei_per_unit"]:
    raise SystemExit("minimum-cell threshold differs from the frozen spec")
if audit["primary_units_observed"] != spec["fixture"]["primary_units_observed"]:
    raise SystemExit("observed primary unit count differs from the frozen spec")
if audit["missing_encoded_as_biological_zero"] is not False:
    raise SystemExit("masked units must never be encoded as zero")
if audit["reconciliation"]["corgi_agreement"] != "exact":
    raise SystemExit("Corgi cross-check did not reconcile")

freeze_tree(root, {
    "artifact_class": "donor_lineage_rna_pseudobulk",
    "dataset_id": "gse296875",
    "donors": 39,
    "biological_unit": "donor_by_lineage_pseudobulk",
    "unit_of_inference": "donor",
    "campaign_spec_sha256": "474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58",
    "outcomes_read": False,
    "model_training_activated": False,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE296875 donor-by-lineage pseudobulk: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
