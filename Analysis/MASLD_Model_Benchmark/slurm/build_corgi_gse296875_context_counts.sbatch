#!/usr/bin/env bash
#SBATCH --job-name=masld-corgi-context-counts
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PYTHON="${PROJECT}/.agent/scanpy_env/bin/python"
MEMBERSHIP="${ROOT}/executions/gse296875-fragment-membership-21063829"
LOCK="${PROJECT}/Analysis/SingleCell/continual_integration/reference/gse296875_source_lock_v44.json"
CORGI="${ROOT}/executions/upstream-acquisitions/corgi-d84c000-zenodo18630048-21064463"
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
TARGET="${ROOT}/executions/corgi-gse296875-context-counts-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite Corgi context artifact: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
"${PYTHON}" -m unittest -v tests.unit.test_build_corgi_gse296875_context_counts
"${PYTHON}" scripts/build_corgi_gse296875_context_counts.py \
  --membership-root "${MEMBERSHIP}" \
  --source-lock "${LOCK}" \
  --corgi-root "${CORGI}" \
  --gtf "${GTF}" \
  --output "${STAGE}"
"${PYTHON}" -m pip freeze > "${STAGE}/pip-freeze.txt"
sha256sum scripts/build_corgi_gse296875_context_counts.py tests/unit/test_build_corgi_gse296875_context_counts.py slurm/build_corgi_gse296875_context_counts.sbatch > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH=src python - "${STAGE}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root = Path(sys.argv[1]).resolve(strict=True)
audit = json.loads((root / "audit.json").read_text())
if audit["status"] != "pass_raw_context_counts_mapper_fit_blocked" or audit["outcomes_read"]:
    raise SystemExit("Corgi raw context receipt differs")
freeze_tree(root, {"artifact_class": "masked_donor_lineage_rna_context_counts", "model_id": "corgi_regular", "dataset_id": "gse296875", "donors": 39, "biological_unit": "donor_by_lineage_pseudobulk", "outcomes_read": False, "model_training_activated": False, "status": "passed"})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Corgi GSE296875 context counts: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
