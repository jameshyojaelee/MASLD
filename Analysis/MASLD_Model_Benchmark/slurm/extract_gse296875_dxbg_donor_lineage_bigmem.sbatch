#!/usr/bin/env bash
#SBATCH --job-name=model-data-007
#SBATCH --partition=bigmem
#SBATCH --account=nslab
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

# bigmem silently redirects anything under 500G to cpu, so --mem=500G and
# --qos=interactive are both load bearing here, not padding.  interactive caps
# concurrency at four jobs, which is fine for one.

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
CANDIDATES="${PROJECT}/Analysis/SingleCell/continual_integration/candidates"
RDS="${CANDIDATES}/external-gse296875-download-v44-2026-08-18T000000Z/GSE296875_processed.rds"
RDS_SHA256=6af74b7dbdd3839841f90cb7f3e47a3184b4dd3dd02242f6f3bdbd14533cfcab
MEMBERSHIP="${ROOT}/executions/gse296875-fragment-membership-21063829/cell_membership.tsv.gz"
FOLDS="${ROOT}/executions/gse296875-donor-folds-20260825/split_lock/donor_folds.tsv"
PSEUDOBULK="${ROOT}/executions/gse296875-donor-lineage-pseudobulk-21107249"
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
TARGET="${ROOT}/executions/gse296875-dxbg-donor-lineage-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite dXbg artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
mkdir -p "${STAGE}/validation"

cat > "${STAGE}/validation/bound_fixtures.sha256" <<'SUMS'
582387b78ec25215d45820c24aace12f130b826e91755bf19ac4581730611e64  executions/gse296875-fragment-membership-21063829/ARTIFACTS.json
e6bc9161acaab508c298c49036af1ad10a720949178a14cba42282933195d416  executions/gse296875-donor-folds-20260825/ARTIFACTS.json
474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58  executions/gse296875-phenotype-campaign-spec-21107179/ARTIFACTS.json
e381f60328645d8e0f2887976021094be35f3d5f8d703a78efe9f5ea7fadb1b7  executions/gse296875-donor-lineage-pseudobulk-21107249/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

module load R/4.4-custom
Rscript scripts/extract_gse296875_dxbg_donor_lineage.R \
  --rds "${RDS}" \
  --rds-sha256 "${RDS_SHA256}" \
  --membership "${MEMBERSHIP}" \
  --folds "${FOLDS}" \
  --raw-features "${PSEUDOBULK}/pseudobulk/genes.tsv" \
  --output "${STAGE}/dxbg" \
  > "${STAGE}/validation/extract_stdout.json" 2> "${STAGE}/validation/extract_stderr.txt"

# Convert to the same on-disk layout as the raw fixture so the fitting and
# scoring jobs bind it with no code change.  A structurally missing value stays
# NaN across the conversion and never becomes a zero.
"${PYTHON}" - "${STAGE}/dxbg" <<'PY'
import csv
import sys
from pathlib import Path

import numpy as np

root = Path(sys.argv[1])
with open(root / "units.tsv", encoding="utf-8", newline="") as handle:
    units = list(csv.DictReader(handle, delimiter="\t"))

matrix = np.loadtxt(
    root / "dxbg_counts.tsv.gz", delimiter="\t", skiprows=1,
    usecols=range(1, 36_602), dtype=np.float64,
)
if matrix.shape != (len(units), 36_601):
    raise SystemExit(f"dXbg matrix shape differs: {matrix.shape}")
if (matrix == 0.0).all(axis=1).any():
    raise SystemExit("an all-zero unit row would misencode structural missingness")

np.savez_compressed(
    root / "donor_lineage_counts.npz",
    counts=matrix,
    unit_observed=np.asarray(
        [row["rna_state"] == "observed" for row in units], dtype=np.bool_
    ),
    nuclei=np.asarray([int(row["nuclei"]) for row in units], dtype=np.int64),
)
print(
    f"converted {matrix.shape[0]} units x {matrix.shape[1]} features; "
    f"structurally missing cells: {int(np.isnan(matrix).sum())}"
)
PY

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true

sha256sum \
  scripts/extract_gse296875_dxbg_donor_lineage.R \
  slurm/extract_gse296875_dxbg_donor_lineage_bigmem.sbatch \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "dxbg" / "audit.json").read_text())
if audit["outcomes_read"] or audit["model_training_activated"]:
    raise SystemExit("the dXbg extraction must read no outcome")
if audit["unmapped_are_missing_never_zero"] is not True:
    raise SystemExit("unmapped features must be missing, never zero")
if audit["cell_axis_matched_frozen_roster"] is not True:
    raise SystemExit("the cell axis did not match the frozen roster")
if audit["dxbg_mapped_fraction"] < 0.90:
    raise SystemExit("the feature crosswalk is not clean")
freeze_tree(root, {
    "artifact_class": "donor_lineage_dxbg_decontaminated_counts",
    "dataset_id": "gse296875",
    "donors": 39,
    "assay": "dXbg",
    "decontaminated": True,
    "campaign_spec_sha256": "474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58",
    "outcomes_read": False,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE296875 dXbg donor-by-lineage counts: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
