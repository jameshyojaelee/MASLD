#!/bin/bash
# WGCNA Preparation Script (Step 00)
# Uses the 'masld_unified' environment

MAMBA_EXE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/in-house_MCD_RNAseq/bin/bin/micromamba"
ENV_PATH="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/unified_env"
SCRIPT="00_prepare_multidataset.R"

echo "--- Starting WGCNA Step 00: Data Preparation ---"
echo "Using Environment: $ENV_PATH"

$MAMBA_EXE run -p $ENV_PATH Rscript $SCRIPT

echo "--- Data Preparation (Step 00) Complete ---"
