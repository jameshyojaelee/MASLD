#!/bin/bash
# run_sex_v2_downstream.sh
# Orchestrator: re-run all downstream scripts affected by the sex v2
# interaction-based classification change.
#
# Usage: bash scripts/run_sex_v2_downstream.sh

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
LOGDIR="${BASE}/RNA-seq/results/stratified_causal/logs"
TMPDIR="${BASE}/RNA-seq/results/stratified_causal/logs/sbatch_tmp"
mkdir -p "$LOGDIR" "$TMPDIR"

echo "============================================================"
echo "  SEX V2 DOWNSTREAM PIPELINE"
echo "  Date: $(date)"
echo "============================================================"

# Helper: write a temporary sbatch script and submit it
submit_script() {
  local name="$1"
  local rscript="$2"
  local dep="${3:-}"
  local script_file="${TMPDIR}/${name}.sbatch"

  cat > "$script_file" << SBEOF
#!/bin/bash
#SBATCH --job-name=${name}
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=${LOGDIR}/${name}_%j.out
#SBATCH --error=${LOGDIR}/${name}_%j.err

set -euo pipefail
eval "\$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "=== ${name}: \$(date) ==="
cd ${BASE}/RNA-seq
Rscript ${rscript}
echo "=== ${name} DONE: \$(date) ==="
SBEOF

  if [ -n "$dep" ]; then
    sbatch --parsable --dependency="afterok:${dep}" "$script_file"
  else
    sbatch --parsable "$script_file"
  fi
}

# ---- STAGE 1: Atlas assembly ----
echo ""
echo "=== STAGE 1: Atlas assembly (27a) ==="
JOB1=$(submit_script "sex_27a" "27a_assemble_evidence_atlas.R")
echo "  27a submitted: ${JOB1}"

# ---- STAGE 2: Parallel analysis (after 27a) ----
echo ""
echo "=== STAGE 2: Parallel analysis scripts (after 27a) ==="

JOB_207=$(submit_script "sex_207" "207_sex_stratified_coloc.R" "${JOB1}")
echo "  207 submitted: ${JOB_207} (dep: ${JOB1})"

JOB_208=$(submit_script "sex_208" "208_subtype_coloc.R" "${JOB1}")
echo "  208 submitted: ${JOB_208} (dep: ${JOB1})"

JOB_215=$(submit_script "sex_215" "215_pharmacogenomic_targets.R" "${JOB1}")
echo "  215 submitted: ${JOB_215} (dep: ${JOB1})"

JOB_202=$(submit_script "sex_202" "202_geneset_enrichment.R" "${JOB1}")
echo "  202 submitted: ${JOB_202} (dep: ${JOB1})"

# ---- STAGE 3: Figures (after Stage 2) ----
echo ""
echo "=== STAGE 3: Figure scripts ==="

JOB_209=$(submit_script "sex_209" "209_sex_subtype_figures.R" "${JOB_207}:${JOB_208}")
echo "  209 submitted: ${JOB_209} (dep: ${JOB_207},${JOB_208})"

JOB_216=$(submit_script "sex_216" "216_pharma_figures.R" "${JOB_215}")
echo "  216 submitted: ${JOB_216} (dep: ${JOB_215})"

# ---- STAGE 4: Final atlas integration (after all Stage 2) ----
echo ""
echo "=== STAGE 4: Final atlas integration (217) ==="

JOB_217=$(submit_script "sex_217" "217_stratified_causal_atlas.R" "${JOB_207}:${JOB_208}:${JOB_215}")
echo "  217 submitted: ${JOB_217} (dep: ${JOB_207},${JOB_208},${JOB_215})"

echo ""
echo "============================================================"
echo "  ALL JOBS SUBMITTED"
echo "  Stage 1: ${JOB1} (27a atlas)"
echo "  Stage 2: ${JOB_207} (207), ${JOB_208} (208), ${JOB_215} (215), ${JOB_202} (202)"
echo "  Stage 3: ${JOB_209} (209 figs), ${JOB_216} (216 figs)"
echo "  Stage 4: ${JOB_217} (217 final atlas)"
echo "============================================================"
echo ""
echo "Monitor: squeue -u \$USER --name=sex_27a,sex_207,sex_208,sex_215,sex_202,sex_209,sex_216,sex_217"
