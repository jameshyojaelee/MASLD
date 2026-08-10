#!/usr/bin/env bash
#SBATCH --job-name=chromBPNet
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/promote_currin_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2/logs/promote_currin_%j.err
set -euo pipefail

# Phase 0 (2026-07-16): promote the AUDITED-LEGITIMATE Currin orientation repair into
# canonical. The original truth set beta_alt=|beta| (constant-sign -> at-chance direction
# by construction). The fixed calibration/01_build_currin_truth.py retains the signed
# FastQTL ALT-dosage beta (model-blind, EA/NEA cross-checked). This regenerates canonical
# truth + calibration WITH THE FIXED SCRIPT (not by copying the isolated repair tree),
# reusing the UNTOUCHED model scores, then finalizes the ChromBPNet gate. It asserts the
# 5 canonical APPLY-ONLY firewall files are byte-unchanged.

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
ENV=${ROOT}/.mamba/seqfunc_chrombpnet
SRC=${ROOT}/GWAS/finemapping/src/seqfunc_v2
BASE=${ROOT}/GWAS/finemapping/results/seqfunc/adult_liver_chrombpnet/hepatocyte_5fold_v2
CANONICAL_TRUTH=${ROOT}/GWAS/finemapping/results/seqfunc/chrombpnet_caqtl/v2_truth
CAL=${BASE}/currin_calibration
FIREWALL=${ROOT}/GWAS/finemapping/results/seqfunc/releases/2026-07-13-r1/canonical_firewall_before.tsv
BACKUP=${BASE}/currin_superseded_preorientation_20260716
PROV=${BASE}/currin_promotion_20260716

mkdir -p "${PROV}" "${BASE}/logs"
test ! -e "${BACKUP}"   # never clobber an existing backup

# --- 1. Back up the canonical files this job overwrites (reversibility + provenance) ---
mkdir -p "${BACKUP}/v2_truth"
cp -p "${CANONICAL_TRUTH}"/*.tsv.gz "${CANONICAL_TRUTH}/truth_contract.json" "${BACKUP}/v2_truth/"
cp -p "${CAL}/currin_external_calibration.json" "${BACKUP}/currin_external_calibration.json"
cp -p "${BASE}/gate_verdict.json" "${BACKUP}/gate_verdict.json"

# --- 2. Firewall snapshot BEFORE (recompute the 5 protected hashes) ---
awk -F '\t' 'NR>1{print $1}' "${FIREWALL}" | while read -r p; do sha256sum "$p"; done \
  > "${PROV}/firewall_before.sha256"

# --- 3. Rebuild canonical truth with the FIXED script (signed ALT-dosage beta) ---
"${ENV}/bin/python" "${SRC}/calibration/01_build_currin_truth.py" \
  --root "${ROOT}" --out "${CANONICAL_TRUTH}"

# --- 4. Regenerate canonical calibration inputs from the corrected truth ---
"${ENV}/bin/python" "${SRC}/chrombpnet/10_prepare_currin.py" \
  --base "${BASE}" --truth "${CANONICAL_TRUTH}" --calibration-dir "${CAL}"

# --- 5. Re-calibrate against the UNTOUCHED model scores (scores dir not rewritten) ---
"${ENV}/bin/python" "${SRC}/chrombpnet/12_calibrate_currin.py" \
  --base "${BASE}" --calibration-dir "${CAL}" \
  --scores-dir "${CAL}/scores" \
  --peak-file "${BASE}/hepatocyte.idr_pooled_summit.narrowPeak" \
  --bootstraps 1000 --workers "${SLURM_CPUS_PER_TASK}"

# --- 6. Finalize the canonical ChromBPNet gate from internal QC + corrected calibration ---
"${ENV}/bin/python" "${SRC}/chrombpnet/06_update_gate.py" \
  --internal "${BASE}/fold0/internal_qc.json" \
  --external "${CAL}/currin_external_calibration.json"

# --- 7. Firewall snapshot AFTER; HARD-assert the 5 protected files are byte-unchanged ---
awk -F '\t' 'NR>1{print $1}' "${FIREWALL}" | while read -r p; do sha256sum "$p"; done \
  > "${PROV}/firewall_after.sha256"
cmp "${PROV}/firewall_before.sha256" "${PROV}/firewall_after.sha256" \
  && echo "FIREWALL OK: 5/5 canonical files byte-unchanged"

# --- 8. Verify: corrected truth counts preserved + gate now unlocked ---
"${ENV}/bin/python" - "${CANONICAL_TRUTH}/truth_contract.json" "${BASE}/gate_verdict.json" <<'PY'
import json, sys
truth = json.load(open(sys.argv[1])); gate = json.load(open(sys.argv[2]))
assert truth["n_positive_1kb"] == 35361, truth["n_positive_1kb"]
assert truth["n_source_background"] == 125689, truth["n_source_background"]
assert truth["one_mb_signed_use"] is False
pos = truth.get("n_beta_alt_positive_1kb"); neg = truth.get("n_beta_alt_negative_1kb")
assert pos and neg and pos > 0 and neg > 0, (pos, neg)   # both signs present (bug guard)
assert gate["currin_external_calibration_pass"] is True, gate
assert gate["saturation_unlocked"] is True, gate
assert gate["status"] == "pass", gate
print("TRUTH counts preserved:", truth["n_positive_1kb"], "pos/neg beta_alt:", pos, neg)
print("GATE:", json.dumps({k: gate[k] for k in ("currin_external_calibration_pass","saturation_unlocked","status")}))
print("PROMOTION COMPLETE")
PY
