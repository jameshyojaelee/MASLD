#!/usr/bin/env bash
#SBATCH --job-name=masld-bulk-activation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

set -euo pipefail
umask 027
export PYTHONDONTWRITEBYTECODE=1
export PYTHONNOUSERSITE=1
: "${COHORT:?COHORT must be exported}"

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SOURCE="${ROOT}/executions/bulk-expansion-source-audit-21065762"
MATRIX="${ROOT}/executions/bulk-expansion-matrix-audit-21065853"
GTF=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
BASE="${PROJECT}/.agent/scanpy_env/bin/python"
TARGET="${ROOT}/executions/${COHORT,,}-activation-${SLURM_JOB_ID}"
STAGE="${TARGET}.staging"
FAILED="${TARGET}.failed"
preserve_failed_stage() {
  if [[ -d "${STAGE}" && ! -e "${FAILED}" ]]; then mv "${STAGE}" "${FAILED}"; fi
}
trap preserve_failed_stage EXIT
if [[ -e "${TARGET}" || -e "${STAGE}" || -e "${FAILED}" ]]; then
  printf 'Refusing to overwrite bulk activation: %s\n' "${TARGET}" >&2
  exit 1
fi
cd "${ROOT}"
"${BASE}" -m unittest -v tests.unit.test_audit_bulk_expansion_activation
"${BASE}" scripts/audit_bulk_expansion_activation.py \
  --cohort "${COHORT}" \
  --source "${SOURCE}" --source-artifacts-sha256 883f870e2a76f78c39d257c24f8e795d2e74f1bd32aba785251b52aa7cb8229e \
  --matrix "${MATRIX}" --matrix-artifacts-sha256 3f934ab9c6ee730c6bd0887b10d963f8210f274c36dc0878cb814c85026dc038 \
  --gtf "${GTF}" --gtf-sha256 73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9 \
  --output "${STAGE}"
sha256sum scripts/audit_bulk_expansion_activation.py tests/unit/test_audit_bulk_expansion_activation.py slurm/audit_bulk_expansion_activation.sbatch > "${STAGE}/source.sha256"
module load python/3.11.5-GCCcore-13.2.0
PYTHONPATH="src" python - "${STAGE}" "${COHORT}" <<'PY'
from pathlib import Path
import json
import sys
from masld_bench.artifacts import freeze_tree

root, cohort = Path(sys.argv[1]).resolve(strict=True), sys.argv[2]
audit = json.loads((root / "activation_audit.json").read_text())
cohort_family_ids = {
    "GSE260666": "gse260666_bulk_rna",
    "GSE268273": "gse268273_imid_masld",
    "GSE274114": "gse274114_mash_hbv",
}
freeze_tree(root, {"artifact_class": "bulk_expansion_activation_audit", "cohort_family_ids": [cohort_family_ids[cohort]], "biological_unit": audit["biological_unit"], "model_training_activated": False})
PY
trap - EXIT
mv "${STAGE}" "${TARGET}"
printf 'Bulk activation audit: %s\n' "${TARGET}"
sha256sum "${TARGET}/ARTIFACTS.json"
