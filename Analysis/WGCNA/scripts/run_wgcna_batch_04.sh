#!/bin/bash
# WGCNA Batch Run Script (Step 04 - Module-Trait Relationship)
# Uses the 'masld_unified' environment

MAMBA_EXE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/in-house_MCD_RNAseq/bin/bin/micromamba"
ENV_PATH="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/unified_env"
SCRIPT="04_module_trait_relationship.R"

echo "--- Starting WGCNA Step 04: Module-Trait Association ---"
echo "Using Environment: $ENV_PATH"

run_step() {
    local dataset=$1
    echo "Running Module-Trait Analysis for $dataset..."
    $MAMBA_EXE run -p $ENV_PATH Rscript $SCRIPT "$dataset"
}

run_step "Patient"
run_step "Other_MCD"
run_step "InHouse_MCD"

echo "--- Step 04 Complete ---"
