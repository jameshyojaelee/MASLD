#!/usr/bin/env bash
##############################################################################
# run_spatial_pipeline.sh — SLURM orchestrator for spatial transcriptomics
#
# Chains all analysis phases with --dependency=afterok.
# Phase 0 (data download, env setup) must be run manually first.
#
# Environment: All scripts use `micromamba activate spatial`
# (NOT rapids_singlecell — spatial env has squidpy, cell2location, etc.)
#
# Usage:
#   cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#   bash Analysis/Spatial/scripts/run_spatial_pipeline.sh
##############################################################################

set -euo pipefail

# ── Parse Arguments ────────────────────────────────────────────────────────
# Usage: bash run_spatial_pipeline.sh [DATASET]
#   DATASET: GSE192741 (default), HRA007511_HMSMA, or ALL
DATASET="${1:-GSE192741}"
SUFFIX=""
if [[ "${DATASET}" == "HRA007511_HMSMA" ]]; then
    SUFFIX="_hmsma"
elif [[ "${DATASET}" == "ALL" ]]; then
    SUFFIX="_all"
fi

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS="Analysis/Spatial/scripts"
LOGS="Analysis/Spatial/logs"
mkdir -p "${LOGS}"

echo "=============================================="
echo "  Spatial Transcriptomics Pipeline"
echo "  Dataset: ${DATASET}"
echo "  Date: $(date)"
echo "=============================================="

# ── Phase 0: Prerequisites Check ────────────────────────────────────────────
echo ""
echo "=== Phase 0: Prerequisites ==="

# Check sample lists exist
for DS in GSE192741 HRA007511 HRA007511_HMSMA; do
    LIST="Analysis/Spatial/metadata/${DS}_samples.txt"
    if [ -f "${LIST}" ]; then
        N=$(wc -l < "${LIST}")
        echo "  ${DS}: ${N} samples"
    fi
done

# Determine which checkpoint to check
CHECKPOINT="Analysis/Spatial/results/preprocessed/merged_spatial${SUFFIX}.h5ad"
if [ ! -f "${CHECKPOINT}" ]; then
    echo ""
    echo "  Preprocessed data not found: ${CHECKPOINT}"
    echo "  Will submit Phase 2 (preprocessing) first."

    # Submit Phase 2
    if [[ "${DATASET}" == "ALL" ]]; then
        JOB_PREPROC=$(sbatch --parsable "${SCRIPTS}/02_build_anndata.sh" \
            --force --output-suffix "${SUFFIX}")
    else
        JOB_PREPROC=$(sbatch --parsable "${SCRIPTS}/02_build_anndata.sh" \
            --force --dataset "${DATASET}" --output-suffix "${SUFFIX}")
    fi
    echo "  Preprocessing: Job ${JOB_PREPROC}"
    PREPROC_DEP="--dependency=afterok:${JOB_PREPROC}"
else
    echo "  Preprocessed data: OK (${CHECKPOINT})"
    PREPROC_DEP=""
fi

# ── Phase 2: Already complete (manual run) ────────────────────────────────
# GSE192741 uses pre-processed SpaceRanger outputs from GEO.
# Phase 2 (02_build_anndata.sh) should be run manually before this orchestrator.

# ── Phase 3: cell2location ──────────────────────────────────────────────────
echo ""
echo "=== Phase 3: cell2location ==="

JOB_REF=$(sbatch --parsable "${SCRIPTS}/03a_prepare_c2l_reference.sh")
echo "  Reference model: Job ${JOB_REF}"

JOB_C2L=$(sbatch --parsable --dependency=afterok:${JOB_REF} \
    "${SCRIPTS}/03b_run_cell2location.sh")
echo "  Spatial model: Job ${JOB_C2L}"

JOB_VAL=$(sbatch --parsable --dependency=afterok:${JOB_C2L} \
    "${SCRIPTS}/03c_validate_deconv.sh")
echo "  Validation: Job ${JOB_VAL}"

# ── Phase 4-5: Analysis (parallel where possible) ───────────────────────────
echo ""
echo "=== Phase 4-5: Analysis ==="

# 4a: Zonation (depends on cell2location)
JOB_ZON=$(sbatch --parsable --dependency=afterok:${JOB_C2L} \
    "${SCRIPTS}/04a_define_zonation.sh")
echo "  Zonation: Job ${JOB_ZON}"

# 4b: DEG zonation mapping (depends on zonation)
JOB_ZON_MAP=$(sbatch --parsable --dependency=afterok:${JOB_ZON} \
    "${SCRIPTS}/04b_deg_zonation_mapping.sh")
echo "  DEG Zonation Map: Job ${JOB_ZON_MAP}"

# 5a: Neighborhood enrichment (depends on cell2location)
JOB_COM=$(sbatch --parsable --dependency=afterok:${JOB_C2L} \
    "${SCRIPTS}/05a_neighborhood_enrichment.sh")
echo "  Neighborhood: Job ${JOB_COM}"

# 5b: Ligand-receptor (depends on cell2location)
JOB_LR=$(sbatch --parsable --dependency=afterok:${JOB_C2L} \
    "${SCRIPTS}/05b_ligand_receptor.sh")
echo "  Ligand-receptor: Job ${JOB_LR}"

# 5c: SVGs (depends only on preprocessing — can run now)
JOB_SVG=$(sbatch --parsable "${SCRIPTS}/05c_spatially_variable_genes.sh")
echo "  SVGs: Job ${JOB_SVG}"

# 5d: Spatial domains (depends on cell2location)
JOB_DOM=$(sbatch --parsable --dependency=afterok:${JOB_C2L} \
    "${SCRIPTS}/05d_spatial_domains.sh")
echo "  Domains: Job ${JOB_DOM}"

# 5e: Co-expression (depends only on preprocessing — can run now)
JOB_COE=$(sbatch --parsable "${SCRIPTS}/05e_spatial_coexpression.sh")
echo "  Co-expression: Job ${JOB_COE}"

# 5f: Trajectory (depends on cell2location)
JOB_TRJ=$(sbatch --parsable --dependency=afterok:${JOB_C2L} \
    "${SCRIPTS}/05f_spatial_trajectory.sh")
echo "  Trajectory: Job ${JOB_TRJ}"

# ── Phase 6: Integration ────────────────────────────────────────────────────
echo ""
echo "=== Phase 6: Integration ==="

JOB_INT=$(sbatch --parsable \
    --dependency=afterok:${JOB_ZON_MAP}:${JOB_COM}:${JOB_LR}:${JOB_SVG}:${JOB_DOM}:${JOB_COE}:${JOB_TRJ}:${JOB_VAL} \
    "${SCRIPTS}/06_integration.sh")
echo "  Integration: Job ${JOB_INT}"

# ── Phase 7: Figures ──────────────────────────────────────────────────────────
echo ""
echo "=== Phase 7: Figures ==="

JOB_FIG=$(sbatch --parsable --dependency=afterok:${JOB_INT} \
    "${SCRIPTS}/07_spatial_figures.sh")
echo "  Figures: Job ${JOB_FIG}"

# ── Summary ─────────────────────────────────────────────────────────────────
echo ""
echo "=============================================="
echo "  Pipeline Submitted"
echo ""
echo "  Dependency DAG:"
echo "    c2l_ref → c2l_spatial ─┬→ Zonation → DEG_Map ─┐"
echo "                           ├→ Neighborhood ────────┤"
echo "                           ├→ Ligand-Receptor ─────┤"
echo "                           ├→ Domains ─────────────├→ Integration → Figures"
echo "                           ├→ Trajectory ──────────┤"
echo "                           └→ Validation ──────────┘"
echo "    (preprocessed) ────────┬→ SVGs ────────────────┘"
echo "                           └→ Co-expression ───────┘"
echo ""
echo "  Monitor: squeue -u \$USER"
echo "  Logs: ${LOGS}/"
echo "=============================================="
