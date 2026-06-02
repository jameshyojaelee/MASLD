#!/bin/bash
# =============================================================================
# run_network_remediation.sh
# -----------------------------------------------------------------------------
# Remediation chain after the gene-ID mismatch bug was identified in the
# multi-evidence network edges. Runs 273 (ID normalization) first, then
# re-runs the downstream pipeline from 262c through 270 using the existing
# sbatch files in logs/network_pipeline/scripts/.
#
# Chain:
#   273_idnorm -> 262c_confidence -> 263_noisyor -> 264_diagnost ->
#   265_assemble -> 266_leiden -> 267_neighbor -> 268_enrich ->
#   269_export -> 270_atlas
#
# All dependencies are `afterok`. No jobs are launched by sourcing this file —
# it submits via `sbatch --parsable`. Run it when ready.
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
LOGDIR="${BASE}/logs/network_pipeline"
SCRIPTDIR="${LOGDIR}/scripts"
mkdir -p "${LOGDIR}" "${SCRIPTDIR}"

echo "=== Network Remediation Chain ==="
echo "Start: $(date)"
echo "Logs:  ${LOGDIR}/"
echo ""

# -----------------------------------------------------------------------------
# Step 0: write the 273 sbatch (doesn't exist in the existing pipeline).
# -----------------------------------------------------------------------------
cat > "${SCRIPTDIR}/273_idnorm.sbatch" << 'EOF'
#!/bin/bash
#SBATCH --job-name=net_273_idnorm
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "273 ID normalization | Start: $(date)"
python 273_normalize_ids.py
echo "273 ID normalization | Exit: $? | End: $(date)"
EOF

# -----------------------------------------------------------------------------
# Helper: submit an existing sbatch file with an optional afterok dependency.
# Usage: submit_chain <step_id> <sbatch_filename> <dep_job_id_or_empty>
# -----------------------------------------------------------------------------
submit_chain() {
  local STEP="$1"
  local SBATCH_FILE="${SCRIPTDIR}/$2"
  local DEP="${3:-}"

  if [[ ! -f "${SBATCH_FILE}" ]]; then
    echo "ERROR: sbatch file not found: ${SBATCH_FILE}" >&2
    exit 1
  fi

  local DEP_FLAG=""
  [[ -n "${DEP}" ]] && DEP_FLAG="--dependency=afterok:${DEP}"

  sbatch --parsable ${DEP_FLAG} \
    --output="${LOGDIR}/${STEP}_%j.out" \
    --error="${LOGDIR}/${STEP}_%j.err" \
    "${SBATCH_FILE}"
}

# -----------------------------------------------------------------------------
# Submit chain
# -----------------------------------------------------------------------------
J273=$(submit_chain 273_idnorm   273_idnorm.sbatch      "")
echo "273 idnorm:       ${J273}  (no dep)"

J262c=$(submit_chain 262c_conf   262c_confidence.sbatch "${J273}")
echo "262c confidence:  ${J262c}  (after 273)"

J263=$(submit_chain 263_noisyor  263_noisyor.sbatch     "${J262c}")
echo "263 noisy-or:     ${J263}  (after 262c)"

J264=$(submit_chain 264_diagnost 264_diagnost.sbatch    "${J263}")
echo "264 diagnostics:  ${J264}  (after 263)"

J265=$(submit_chain 265_assemble 265_assemble.sbatch    "${J264}")
echo "265 assemble:     ${J265}  (after 264)"

J266=$(submit_chain 266_leiden   266_leiden.sbatch      "${J265}")
echo "266 leiden:       ${J266}  (after 265)"

J267=$(submit_chain 267_neighbor 267_neighbor.sbatch    "${J266}")
echo "267 neighborhoods:${J267}  (after 266)"

J268=$(submit_chain 268_enrich   268_enrich.sbatch      "${J266}")
echo "268 enrichment:   ${J268}  (after 266)"

J269=$(submit_chain 269_export   269_export.sbatch      "${J267}:${J268}")
echo "269 export:       ${J269}  (after 267+268)"

J270=$(submit_chain 270_atlas    270_atlas.sbatch       "${J268}")
echo "270 atlas:        ${J270}  (after 268)"

echo ""
echo "=== Remediation chain submitted (10 jobs) ==="
echo "Monitor: squeue -u \$USER | grep net_"
echo "Logs:    ${LOGDIR}/"
echo "End:     $(date)"
