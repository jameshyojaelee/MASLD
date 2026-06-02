#!/bin/bash
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
ENV_PREFIX="/gpfs/commons/groups/sanjana_lab/Cas13/Pathway/.mamba/pathway_analysis"
export LD_LIBRARY_PATH="${ENV_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
$MICROMAMBA run -p "${ENV_PREFIX}" Rscript --vanilla -e ".libPaths(.libPaths()[!grepl('jameslee', .libPaths())])" "$@"
