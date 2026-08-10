#!/usr/bin/env bash
#SBATCH --job-name=bootstrap
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/currin_orientation_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/currin_orientation_%j.err
set -euo pipefail

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
SRC=${ROOT}/GWAS/finemapping/src/seqfunc_v2
BASE=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2
CANONICAL_TRUTH=${ROOT}/GWAS/finemapping/results/seqfunc/chrombpnet_caqtl/v2_truth
REPAIR_ROOT=${BASE}/currin_orientation_repair_20260715
REPAIR_TRUTH=${REPAIR_ROOT}/truth
REPAIR_CAL=${REPAIR_ROOT}/calibration

test ! -e "${REPAIR_ROOT}"
mkdir -p "${REPAIR_ROOT}"

sha256sum \
  "${CANONICAL_TRUTH}/currin_positive_1kb.tsv.gz" \
  "${CANONICAL_TRUTH}/truth_contract.json" \
  "${BASE}/currin_calibration/currin_external_calibration.json" \
  "${BASE}/gate_verdict.json" > "${REPAIR_ROOT}/canonical_before.sha256"

"${ENV}/bin/python" "${SRC}/calibration/01_build_currin_truth.py" \
  --root "${ROOT}" --out "${REPAIR_TRUTH}"

"${ENV}/bin/python" "${SRC}/chrombpnet/10_prepare_currin.py" \
  --base "${BASE}" --truth "${REPAIR_TRUTH}" --calibration-dir "${REPAIR_CAL}"

"${ENV}/bin/python" "${SRC}/chrombpnet/12_calibrate_currin.py" \
  --base "${BASE}" --calibration-dir "${REPAIR_CAL}" \
  --scores-dir "${BASE}/currin_calibration/scores" \
  --peak-file "${BASE}/hepatocyte.idr_pooled_summit.narrowPeak" \
  --bootstraps 1000 --workers "${SLURM_CPUS_PER_TASK}"

sha256sum \
  "${CANONICAL_TRUTH}/currin_positive_1kb.tsv.gz" \
  "${CANONICAL_TRUTH}/truth_contract.json" \
  "${BASE}/currin_calibration/currin_external_calibration.json" \
  "${BASE}/gate_verdict.json" > "${REPAIR_ROOT}/canonical_after.sha256"
cmp "${REPAIR_ROOT}/canonical_before.sha256" "${REPAIR_ROOT}/canonical_after.sha256"
test -s "${REPAIR_CAL}/currin_external_calibration.json"
