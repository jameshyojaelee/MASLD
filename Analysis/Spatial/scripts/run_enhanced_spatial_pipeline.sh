#!/usr/bin/env bash
###############################################################################
# run_enhanced_spatial_pipeline.sh — Master orchestrator for 4 spatial methods
#
# Submits all Phase 0–6 jobs with SLURM dependency chaining.
# Run from login node: bash Analysis/Spatial/scripts/run_enhanced_spatial_pipeline.sh
###############################################################################
set -euo pipefail

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Spatial/scripts
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

echo "============================================================"
echo "  Enhanced Spatial Pipeline — Master Orchestrator"
echo "============================================================"

# Ensure output directories exist
mkdir -p Analysis/Spatial/logs
mkdir -p Analysis/Spatial/data/gsmap_gwas
mkdir -p Analysis/Spatial/data/gsmap_input
mkdir -p Analysis/Spatial/results/{multisp,ontrac,ontrac/input,gsmap,commot,commot/validation_vu,spatial_integration}
mkdir -p scripts/figures/output

###############################################################################
# Phase 0 — Environment setup (parallel)
###############################################################################
echo ""
echo "=== Phase 0: Environment Setup ==="

# Check if spatial_multiomics env already exists
if micromamba env list 2>/dev/null | grep -q spatial_multiomics; then
    echo "  spatial_multiomics env already exists, skipping setup"
    ENV_MULTIOMICS_JOB=""
else
    ENV_MULTIOMICS_JOB=$(sbatch --parsable "$SCRIPT_DIR/00i_setup_spatial_multiomics_env.sh")
    echo "  spatial_multiomics env: SLURM $ENV_MULTIOMICS_JOB (uses conda clone)"
fi

# Check if gsmap env already exists
if micromamba env list 2>/dev/null | grep -q gsmap; then
    echo "  gsmap env already exists, skipping setup"
    ENV_GSMAP_JOB=""
else
    ENV_GSMAP_JOB=$(sbatch --parsable "$SCRIPT_DIR/15a_setup_gsmap_env.sh")
    echo "  gsmap env: SLURM $ENV_GSMAP_JOB"
fi

# Install COMMOT into existing spatial env
COMMOT_INSTALL_JOB=$(sbatch --parsable "$SCRIPT_DIR/16a_install_commot.sh")
echo "  COMMOT install: SLURM $COMMOT_INSTALL_JOB"

###############################################################################
# Phase 1 — Data preparation (parallel, after Phase 0)
###############################################################################
echo ""
echo "=== Phase 1: Data Preparation ==="

# MultiSP: ATAC projection (depends on spatial_multiomics env)
if [ -n "${ENV_MULTIOMICS_JOB:-}" ]; then
    ATAC_PROJ_JOB=$(sbatch --parsable --dependency=afterok:$ENV_MULTIOMICS_JOB "$SCRIPT_DIR/13a_atac_pseudobulk_projection.sh")
else
    ATAC_PROJ_JOB=$(sbatch --parsable "$SCRIPT_DIR/13a_atac_pseudobulk_projection.sh")
fi
echo "  13a ATAC projection: SLURM $ATAC_PROJ_JOB"

# ONTraC: prepare input (depends on spatial_multiomics env)
if [ -n "${ENV_MULTIOMICS_JOB:-}" ]; then
    ONTRAC_PREP_JOB=$(sbatch --parsable --dependency=afterok:$ENV_MULTIOMICS_JOB "$SCRIPT_DIR/14a_ontrac_prepare_input.sh")
else
    ONTRAC_PREP_JOB=$(sbatch --parsable "$SCRIPT_DIR/14a_ontrac_prepare_input.sh")
fi
echo "  14a ONTraC prep: SLURM $ONTRAC_PREP_JOB"

# gsMap: format GWAS (depends on gsmap env)
if [ -n "${ENV_GSMAP_JOB:-}" ]; then
    GSMAP_FMT_JOB=$(sbatch --parsable --dependency=afterok:$ENV_GSMAP_JOB "$SCRIPT_DIR/15b_format_gwas_for_gsmap.sh")
else
    GSMAP_FMT_JOB=$(sbatch --parsable "$SCRIPT_DIR/15b_format_gwas_for_gsmap.sh")
fi
echo "  15b GWAS format: SLURM $GSMAP_FMT_JOB"

# gsMap: prepare spatial data (depends on gsmap env)
if [ -n "${ENV_GSMAP_JOB:-}" ]; then
    GSMAP_PREP_JOB=$(sbatch --parsable --dependency=afterok:$ENV_GSMAP_JOB "$SCRIPT_DIR/15c_prepare_spatial_for_gsmap.sh")
else
    GSMAP_PREP_JOB=$(sbatch --parsable "$SCRIPT_DIR/15c_prepare_spatial_for_gsmap.sh")
fi
echo "  15c Spatial prep: SLURM $GSMAP_PREP_JOB"

# COMMOT: prepare LR database (depends on COMMOT install)
COMMOT_PREP_JOB=$(sbatch --parsable --dependency=afterok:$COMMOT_INSTALL_JOB "$SCRIPT_DIR/16b_prepare_commot.sh")
echo "  16b COMMOT prep: SLURM $COMMOT_PREP_JOB"

###############################################################################
# Phase 2 — Core computation (parallel across methods)
###############################################################################
echo ""
echo "=== Phase 2: Core Computation ==="

# MultiSP domains (depends on 13a)
MULTISP_DOM_JOB=$(sbatch --parsable --dependency=afterok:$ATAC_PROJ_JOB "$SCRIPT_DIR/13b_multisp_domains.sh")
echo "  13b MultiSP domains: SLURM $MULTISP_DOM_JOB"

# ONTraC GNN run (depends on 14a)
ONTRAC_RUN_JOB=$(sbatch --parsable --dependency=afterok:$ONTRAC_PREP_JOB "$SCRIPT_DIR/14b_ontrac_run.sh")
echo "  14b ONTraC run: SLURM $ONTRAC_RUN_JOB"

# gsMap slice mean (depends on 15c)
GSMAP_SLICE_JOB=$(sbatch --parsable --dependency=afterok:$GSMAP_PREP_JOB "$SCRIPT_DIR/15d_create_slice_mean.sh")
echo "  15d Slice mean: SLURM $GSMAP_SLICE_JOB"

# COMMOT run (depends on 16b)
COMMOT_RUN_JOB=$(sbatch --parsable --dependency=afterok:$COMMOT_PREP_JOB "$SCRIPT_DIR/16c_run_commot.sh")
echo "  16c COMMOT run: SLURM $COMMOT_RUN_JOB"

###############################################################################
# Phase 3 — Heavy computation
###############################################################################
echo ""
echo "=== Phase 3: Heavy Computation ==="

# gsMap array (depends on 15a env + 15b GWAS + 15d slice mean)
GSMAP_DEPS="${GSMAP_FMT_JOB}:${GSMAP_SLICE_JOB}"
if [ -n "${ENV_GSMAP_JOB:-}" ]; then
    GSMAP_DEPS="${ENV_GSMAP_JOB}:${GSMAP_DEPS}"
fi
GSMAP_RUN_JOB=$(sbatch --parsable --dependency=afterok:$GSMAP_DEPS "$SCRIPT_DIR/15e_run_gsmap.sh")
echo "  15e gsMap array: SLURM $GSMAP_RUN_JOB"

# ONTraC disease comparison (depends on 14b)
ONTRAC_DISEASE_JOB=$(sbatch --parsable --dependency=afterok:$ONTRAC_RUN_JOB "$SCRIPT_DIR/14c_ontrac_disease_comparison.sh")
echo "  14c ONTraC disease: SLURM $ONTRAC_DISEASE_JOB"

# MultiSP communication (depends on 13b)
MULTISP_COMM_JOB=$(sbatch --parsable --dependency=afterok:$MULTISP_DOM_JOB "$SCRIPT_DIR/13c_multisp_communication.sh")
echo "  13c MultiSP comm: SLURM $MULTISP_COMM_JOB"

# COMMOT direction (depends on 16c)
COMMOT_DIR_JOB=$(sbatch --parsable --dependency=afterok:$COMMOT_RUN_JOB "$SCRIPT_DIR/16d_commot_direction.sh")
echo "  16d COMMOT direction: SLURM $COMMOT_DIR_JOB"

# COMMOT regulated genes (depends on 16c)
COMMOT_REG_JOB=$(sbatch --parsable --dependency=afterok:$COMMOT_RUN_JOB "$SCRIPT_DIR/16e_commot_regulated_genes.sh")
echo "  16e COMMOT regulated: SLURM $COMMOT_REG_JOB"

###############################################################################
# Phase 4 — Analysis & validation (parallel)
###############################################################################
echo ""
echo "=== Phase 4: Analysis & Validation ==="

# MultiSP validation (depends on 13b, 13c)
MULTISP_VAL_JOB=$(sbatch --parsable --dependency=afterok:$MULTISP_DOM_JOB:$MULTISP_COMM_JOB "$SCRIPT_DIR/13d_multisp_validation.sh")
echo "  13d MultiSP validation: SLURM $MULTISP_VAL_JOB"

# ONTraC zonation overlay (depends on 14b, 14c)
ONTRAC_ZON_JOB=$(sbatch --parsable --dependency=afterok:$ONTRAC_RUN_JOB:$ONTRAC_DISEASE_JOB "$SCRIPT_DIR/14d_ontrac_zonation_overlay.sh")
echo "  14d ONTraC zonation: SLURM $ONTRAC_ZON_JOB"

# gsMap analysis (depends on 15e)
GSMAP_ANALYSIS_JOB=$(sbatch --parsable --dependency=afterok:$GSMAP_RUN_JOB "$SCRIPT_DIR/15f_gsmap_analysis.sh")
echo "  15f gsMap analysis: SLURM $GSMAP_ANALYSIS_JOB"

# COMMOT differential (depends on 16c, 16d)
COMMOT_DIFF_JOB=$(sbatch --parsable --dependency=afterok:$COMMOT_RUN_JOB:$COMMOT_DIR_JOB "$SCRIPT_DIR/16f_commot_differential.sh")
echo "  16f COMMOT diff: SLURM $COMMOT_DIFF_JOB"

# COMMOT convergence (depends on 16e, 16f)
COMMOT_CONV_JOB=$(sbatch --parsable --dependency=afterok:$COMMOT_REG_JOB:$COMMOT_DIFF_JOB "$SCRIPT_DIR/16g_commot_convergence.sh")
echo "  16g COMMOT convergence: SLURM $COMMOT_CONV_JOB"

# COMMOT Vu validation (depends on 16b — can run in parallel with GSE192741)
COMMOT_VU_JOB=$(sbatch --parsable --dependency=afterok:$COMMOT_PREP_JOB "$SCRIPT_DIR/16h_commot_validation.sh")
echo "  16h COMMOT Vu valid: SLURM $COMMOT_VU_JOB"

###############################################################################
# Phase 5 — Integration
###############################################################################
echo ""
echo "=== Phase 5: Per-Method Integration ==="

# ONTraC integration (depends on 14d)
ONTRAC_INT_JOB=$(sbatch --parsable --dependency=afterok:$ONTRAC_ZON_JOB "$SCRIPT_DIR/14e_ontrac_integration.sh")
echo "  14e ONTraC integrate: SLURM $ONTRAC_INT_JOB"

# gsMap integration (depends on 15f)
GSMAP_INT_JOB=$(sbatch --parsable --dependency=afterok:$GSMAP_ANALYSIS_JOB "$SCRIPT_DIR/15g_integrate_gsmap.sh")
echo "  15g gsMap integrate: SLURM $GSMAP_INT_JOB"

# COMMOT integration (depends on 16g, 16h)
COMMOT_INT_JOB=$(sbatch --parsable --dependency=afterok:$COMMOT_CONV_JOB:$COMMOT_VU_JOB "$SCRIPT_DIR/16i_integrate_commot.sh")
echo "  16i COMMOT integrate: SLURM $COMMOT_INT_JOB"

###############################################################################
# Phase 6 — Final unified integration + figures
###############################################################################
echo ""
echo "=== Phase 6: Unified Integration & Figures ==="

# Unified spatial integration (depends on all per-method integrations + MultiSP validation)
UNIFIED_JOB=$(sbatch --parsable --dependency=afterok:$MULTISP_VAL_JOB:$ONTRAC_INT_JOB:$GSMAP_INT_JOB:$COMMOT_INT_JOB \
    "$SCRIPT_DIR/17a_unified_spatial_integration.sh")
echo "  17a Unified integration: SLURM $UNIFIED_JOB"

# Main spatial figure (depends on unified integration)
FIGURE_JOB=$(sbatch --parsable --dependency=afterok:$UNIFIED_JOB "$SCRIPT_DIR/17b_spatial_main_figure.sh")
echo "  17b Main figure: SLURM $FIGURE_JOB"

###############################################################################
# Summary
###############################################################################
echo ""
echo "============================================================"
echo "  All jobs submitted. Dependency chain:"
echo "============================================================"
echo ""
echo "  Phase 0 (env):     ${ENV_MULTIOMICS_JOB:-skip} ${ENV_GSMAP_JOB:-skip} $COMMOT_INSTALL_JOB"
echo "  Phase 1 (prep):    $ATAC_PROJ_JOB $ONTRAC_PREP_JOB $GSMAP_FMT_JOB $GSMAP_PREP_JOB $COMMOT_PREP_JOB"
echo "  Phase 2 (core):    $MULTISP_DOM_JOB $ONTRAC_RUN_JOB $GSMAP_SLICE_JOB $COMMOT_RUN_JOB"
echo "  Phase 3 (heavy):   $GSMAP_RUN_JOB $ONTRAC_DISEASE_JOB $MULTISP_COMM_JOB $COMMOT_DIR_JOB $COMMOT_REG_JOB"
echo "  Phase 4 (analysis): $MULTISP_VAL_JOB $ONTRAC_ZON_JOB $GSMAP_ANALYSIS_JOB $COMMOT_DIFF_JOB $COMMOT_CONV_JOB $COMMOT_VU_JOB"
echo "  Phase 5 (integrate): $ONTRAC_INT_JOB $GSMAP_INT_JOB $COMMOT_INT_JOB"
echo "  Phase 6 (final):   $UNIFIED_JOB $FIGURE_JOB"
echo ""
echo "  Monitor: squeue -u \$USER"
echo "  Final figure: scripts/figures/output/fig_spatial_main.pdf"
echo "  Last job (figure): $FIGURE_JOB"
echo ""

# Save job IDs for tracking
cat > Analysis/Spatial/logs/pipeline_job_ids.txt << JOBEOF
# Enhanced Spatial Pipeline — Job IDs ($(date))
PHASE0_MULTIOMICS=${ENV_MULTIOMICS_JOB:-skip}
PHASE0_GSMAP=${ENV_GSMAP_JOB:-skip}
PHASE0_COMMOT=$COMMOT_INSTALL_JOB
PHASE1_ATAC=$ATAC_PROJ_JOB
PHASE1_ONTRAC=$ONTRAC_PREP_JOB
PHASE1_GWAS_FMT=$GSMAP_FMT_JOB
PHASE1_SPATIAL_PREP=$GSMAP_PREP_JOB
PHASE1_COMMOT_PREP=$COMMOT_PREP_JOB
PHASE2_MULTISP=$MULTISP_DOM_JOB
PHASE2_ONTRAC=$ONTRAC_RUN_JOB
PHASE2_GSMAP_SLICE=$GSMAP_SLICE_JOB
PHASE2_COMMOT=$COMMOT_RUN_JOB
PHASE3_GSMAP_ARRAY=$GSMAP_RUN_JOB
PHASE3_ONTRAC_DISEASE=$ONTRAC_DISEASE_JOB
PHASE3_MULTISP_COMM=$MULTISP_COMM_JOB
PHASE3_COMMOT_DIR=$COMMOT_DIR_JOB
PHASE3_COMMOT_REG=$COMMOT_REG_JOB
PHASE4_MULTISP_VAL=$MULTISP_VAL_JOB
PHASE4_ONTRAC_ZON=$ONTRAC_ZON_JOB
PHASE4_GSMAP_ANALYSIS=$GSMAP_ANALYSIS_JOB
PHASE4_COMMOT_DIFF=$COMMOT_DIFF_JOB
PHASE4_COMMOT_CONV=$COMMOT_CONV_JOB
PHASE4_COMMOT_VU=$COMMOT_VU_JOB
PHASE5_ONTRAC_INT=$ONTRAC_INT_JOB
PHASE5_GSMAP_INT=$GSMAP_INT_JOB
PHASE5_COMMOT_INT=$COMMOT_INT_JOB
PHASE6_UNIFIED=$UNIFIED_JOB
PHASE6_FIGURE=$FIGURE_JOB
JOBEOF

echo "  Job IDs saved to: Analysis/Spatial/logs/pipeline_job_ids.txt"
