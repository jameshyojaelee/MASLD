#!/bin/bash
# run_progression_analysis.sh
# ---------------------------------------------------------------------------
# SLURM orchestrator for the complete progression analysis pipeline.
# Chains Scripts 130-146 with SLURM dependencies across 4 phases:
#
# Phase 1 (independent): 130, 131, 132, 141 — Core DE + ATAC
# Phase 2 (depends on P1): 133, 134, 139, 140, 142, 144
# Phase 3 (depends on P2): 135, 138, 143
# Phase 4 (depends on all): 145, 146
#
# Usage: bash run_progression_analysis.sh
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPTS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
LOGS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs"
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"

mkdir -p "$LOGS"

# Common SLURM preamble
INIT_CMD='eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)" && micromamba activate rnaseq'

echo "=== Progression Analysis Pipeline ==="
echo "Started: $(date)"
echo ""

# ============================================================
# Phase 1: Core DE (independent jobs)
# ============================================================
echo "--- Phase 1: Core Differential Expression ---"

# Script 130: Binary progression contrasts (16 CPUs, 64G, ~12h)
JOB_130=$(sbatch --parsable \
    --job-name=stg130_contrasts \
    --partition=cpu --cpus-per-task=16 --mem=64G --time=48:00:00 \
    --output="$LOGS/130_progression_contrasts_%j.out" \
    --error="$LOGS/130_progression_contrasts_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/130_progression_contrasts_dream.R'")
echo "  Script 130 (binary contrasts): Job $JOB_130"

# Script 131: NAS component ordinal (8 CPUs, 32G)
JOB_131=$(sbatch --parsable \
    --job-name=stg131_ordinal \
    --partition=cpu --cpus-per-task=8 --mem=32G --time=48:00:00 \
    --output="$LOGS/131_nas_ordinal_%j.out" \
    --error="$LOGS/131_nas_ordinal_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/131_nas_component_ordinal_dream.R'")
echo "  Script 131 (NAS ordinal): Job $JOB_131"

# Script 132: Sex-stratified progression (16 CPUs, 64G)
JOB_132=$(sbatch --parsable \
    --job-name=stg132_sex \
    --partition=cpu --cpus-per-task=16 --mem=64G --time=48:00:00 \
    --output="$LOGS/132_sex_stratified_%j.out" \
    --error="$LOGS/132_sex_stratified_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/132_sex_stratified_progression.R'")
echo "  Script 132 (sex-stratified): Job $JOB_132"

# Script 141: ATAC MASL-vs-MASH (8 CPUs, 64G, snapatac2 env)
JOB_141=$(sbatch --parsable \
    --job-name=stg141_atac \
    --partition=cpu --cpus-per-task=8 --mem=64G --time=48:00:00 \
    --output="$LOGS/141_atac_da_%j.out" \
    --error="$LOGS/141_atac_da_%j.err" \
    --wrap="bash -c 'eval \"\$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)\" && micromamba activate snapatac2 && cd ${BASE} && python ${SCRIPTS}/141_atac_masl_vs_mash.py'")
echo "  Script 141 (ATAC DA): Job $JOB_141"

echo ""

# ============================================================
# Phase 2: Downstream re-analysis (depends on Phase 1)
# ============================================================
echo "--- Phase 2: Downstream Re-analysis ---"
P1_DEPS="--dependency=afterok:${JOB_130}:${JOB_131}:${JOB_132}"

# Script 133: Pathway enrichment (8 CPUs, 32G)
JOB_133=$(sbatch --parsable \
    ${P1_DEPS} \
    --job-name=stg133_gsea \
    --partition=cpu --cpus-per-task=8 --mem=32G --time=48:00:00 \
    --output="$LOGS/133_gsea_%j.out" \
    --error="$LOGS/133_gsea_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/133_progression_pathway_enrichment.R'")
echo "  Script 133 (GSEA): Job $JOB_133 (after $JOB_130,$JOB_131,$JOB_132)"

# Script 134: Consensus matrix (4 CPUs, 16G)
JOB_134=$(sbatch --parsable \
    ${P1_DEPS} \
    --job-name=stg134_consensus \
    --partition=cpu --cpus-per-task=4 --mem=32G --time=48:00:00 \
    --output="$LOGS/134_consensus_%j.out" \
    --error="$LOGS/134_consensus_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/134_progression_consensus.R'")
echo "  Script 134 (consensus): Job $JOB_134 (after P1)"

# Script 139: Deconvolution attribution (16 CPUs, 64G)
JOB_139=$(sbatch --parsable \
    --dependency=afterok:${JOB_130} \
    --job-name=stg139_deconv \
    --partition=cpu --cpus-per-task=16 --mem=64G --time=48:00:00 \
    --output="$LOGS/139_deconv_%j.out" \
    --error="$LOGS/139_deconv_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/139_progression_deconv_attribution.R'")
echo "  Script 139 (deconv): Job $JOB_139 (after $JOB_130)"

# Script 140: Cross-species concordance (8 CPUs, 32G)
JOB_140=$(sbatch --parsable \
    --dependency=afterok:${JOB_130} \
    --job-name=stg140_concordance \
    --partition=cpu --cpus-per-task=8 --mem=32G --time=48:00:00 \
    --output="$LOGS/140_concordance_%j.out" \
    --error="$LOGS/140_concordance_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/140_progression_concordance.R'")
echo "  Script 140 (concordance): Job $JOB_140 (after $JOB_130)"

# Script 142: TF activity (8 CPUs, 32G)
JOB_142=$(sbatch --parsable \
    --dependency=afterok:${JOB_130} \
    --job-name=stg142_tf \
    --partition=cpu --cpus-per-task=8 --mem=32G --time=48:00:00 \
    --output="$LOGS/142_tf_activity_%j.out" \
    --error="$LOGS/142_tf_activity_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/142_progression_tf_activity.R'")
echo "  Script 142 (TF activity): Job $JOB_142 (after $JOB_130)"

# Script 144: Subtype × progression (4 CPUs, 16G)
JOB_144=$(sbatch --parsable \
    --dependency=afterok:${JOB_130} \
    --job-name=stg144_subtype \
    --partition=cpu --cpus-per-task=4 --mem=16G --time=48:00:00 \
    --output="$LOGS/144_subtype_%j.out" \
    --error="$LOGS/144_subtype_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/144_subtype_progression_interaction.R'")
echo "  Script 144 (subtype): Job $JOB_144 (after $JOB_130)"

echo ""

# ============================================================
# Phase 3: Causal + drug + regulatory (depends on Phase 2)
# ============================================================
echo "--- Phase 3: Causal Inference & Drug Analysis ---"

# Script 135: Progression TWAS (16 CPUs, 200G, bigmem)
JOB_135=$(sbatch --parsable \
    --dependency=afterok:${JOB_130} \
    --job-name=stg135_twas \
    --partition=bigmem --cpus-per-task=16 --mem=200G --time=48:00:00 \
    --output="$LOGS/135_twas_%j.out" \
    --error="$LOGS/135_twas_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/135_progression_twas.R'")
echo "  Script 135 (TWAS): Job $JOB_135 (bigmem, after $JOB_130)"

# Script 138: Drug repurposing (8 CPUs, 32G)
JOB_138=$(sbatch --parsable \
    --dependency=afterok:${JOB_130} \
    --job-name=stg138_drugs \
    --partition=cpu --cpus-per-task=8 --mem=32G --time=48:00:00 \
    --output="$LOGS/138_drugs_%j.out" \
    --error="$LOGS/138_drugs_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/138_progression_drug_repurposing.R'")
echo "  Script 138 (drugs): Job $JOB_138 (after $JOB_130)"

# Script 143: Regulon convergence (4 CPUs, 16G, depends on 142 + 138)
JOB_143=$(sbatch --parsable \
    --dependency=afterok:${JOB_142}:${JOB_138} \
    --job-name=stg143_convergence \
    --partition=cpu --cpus-per-task=4 --mem=16G --time=48:00:00 \
    --output="$LOGS/143_convergence_%j.out" \
    --error="$LOGS/143_convergence_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/143_progression_regulon_convergence.R'")
echo "  Script 143 (convergence): Job $JOB_143 (after $JOB_142,$JOB_138)"

echo ""

# ============================================================
# Phase 4: Integration & synthesis (depends on all)
# ============================================================
echo "--- Phase 4: Integration & Synthesis ---"
ALL_P2_P3="${JOB_133}:${JOB_134}:${JOB_135}:${JOB_138}:${JOB_139}:${JOB_140}:${JOB_141}:${JOB_142}:${JOB_143}:${JOB_144}"

# Script 145: Atlas expansion (4 CPUs, 32G)
JOB_145=$(sbatch --parsable \
    --dependency=afterok:${ALL_P2_P3} \
    --job-name=stg145_atlas \
    --partition=cpu --cpus-per-task=4 --mem=32G --time=48:00:00 \
    --output="$LOGS/145_atlas_%j.out" \
    --error="$LOGS/145_atlas_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/145_expand_atlas_progression.R'")
echo "  Script 145 (atlas): Job $JOB_145 (after all P2+P3)"

# Script 146: Synthesis figure (4 CPUs, 16G, depends on 145)
JOB_146=$(sbatch --parsable \
    --dependency=afterok:${JOB_145} \
    --job-name=stg146_synthesis \
    --partition=cpu --cpus-per-task=4 --mem=16G --time=48:00:00 \
    --output="$LOGS/146_synthesis_%j.out" \
    --error="$LOGS/146_synthesis_%j.err" \
    --wrap="bash -c '${INIT_CMD} && cd ${BASE} && Rscript ${SCRIPTS}/146_progression_synthesis.R'")
echo "  Script 146 (synthesis): Job $JOB_146 (after $JOB_145)"

echo ""
echo "=== All jobs submitted ==="
echo ""
echo "Phase 1 (independent):  $JOB_130, $JOB_131, $JOB_132, $JOB_141"
echo "Phase 2 (after P1):     $JOB_133, $JOB_134, $JOB_139, $JOB_140, $JOB_142, $JOB_144"
echo "Phase 3 (after P1/P2):  $JOB_135, $JOB_138, $JOB_143"
echo "Phase 4 (after all):    $JOB_145, $JOB_146"
echo ""
echo "Monitor: squeue -u \$USER --format='%.10i %.15j %.8T %.10M %.6C %.10m %.30R'"
echo "Total jobs: 15"
