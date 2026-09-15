#!/usr/bin/env bash
#SBATCH --job-name=model-check-005
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1

: "${PSEUDOBULK:?PSEUDOBULK path is required}"
: "${PSEUDOBULK_SHA256:?PSEUDOBULK_SHA256 is required}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
TARGET="${ROOT}/executions/gse296875-lineage-ambient-load-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"

preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT

if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite ambient audit artifact: %s\n' "${TARGET}" >&2
  exit 1
fi

cd "${ROOT}"
mkdir -p "${STAGE}/validation"

cat > "${STAGE}/validation/bound_fixtures.sha256" <<SUMS
474db3b44b22ff3c3556d1fcca92c5c17d8a10f2b873b5c7ae345c81d6382a58  executions/gse296875-phenotype-campaign-spec-21107179/ARTIFACTS.json
${PSEUDOBULK_SHA256}  ${PSEUDOBULK}/ARTIFACTS.json
SUMS
sha256sum --check --strict "${STAGE}/validation/bound_fixtures.sha256"

"${PYTHON}" scripts/audit_gse296875_lineage_ambient_load.py \
  --pseudobulk "${PSEUDOBULK}" \
  --pseudobulk-sha256 "${PSEUDOBULK_SHA256}" \
  --output "${STAGE}/ambient" \
  > "${STAGE}/validation/ambient_stdout.json"

"${PYTHON}" -m pip freeze --all > "${STAGE}/validation/pip_freeze.txt" 2>&1 || true

sha256sum \
  scripts/audit_gse296875_lineage_ambient_load.py \
  slurm/audit_gse296875_lineage_ambient_load_cpu.sbatch \
  > "${STAGE}/source.sha256"

module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python3 - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "ambient" / "audit.json").read_text())
if audit["outcomes_read"] or audit["model_fitted"]:
    raise SystemExit("the ambient audit must read no outcome and fit no model")
freeze_tree(root, {
    "artifact_class": "donor_lineage_ambient_load_descriptive",
    "dataset_id": "gse296875",
    "outcomes_read": False,
    "model_fitted": False,
    "status": "passed",
})
PY

trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'GSE296875 lineage ambient load: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
