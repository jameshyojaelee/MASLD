#!/bin/bash
# WGCNA Batch Run Script (Step 03)
# Uses the 'masld_unified' environment

MAMBA_EXE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/in-house_MCD_RNAseq/bin/bin/micromamba"
ENV_PATH="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/unified_env"
SCRIPT="03_network_construction.R"
POWER=6

echo "--- Starting WGCNA Batch Run (Power=$POWER) ---"
echo "Using Environment: $ENV_PATH"

run_wgcna() {
    local dataset=$1
    echo "Running $dataset..."
    $MAMBA_EXE run -p $ENV_PATH Rscript $SCRIPT "$dataset" $POWER
}

run_wgcna "Patient"
run_wgcna "Other_MCD"
run_wgcna "InHouse_MCD"

echo "--- Batch Run Complete ---"
