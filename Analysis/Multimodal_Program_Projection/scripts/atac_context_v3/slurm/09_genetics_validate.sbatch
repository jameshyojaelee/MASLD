#!/bin/bash
#SBATCH --job-name=atacgenval
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_genvalidate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Multimodal_Program_Projection/logs/atac_context_v3_genvalidate_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PACKAGE=${ROOT}/Analysis/Multimodal_Program_Projection/scripts/atac_context_v3
MAMBA=/gpfs/commons/home/jameslee/.local/bin/micromamba
eval "$(${MAMBA} shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u
Rscript "${PACKAGE}/08_liftover_variants.R"
set +u
micromamba activate snapatac2
set -u
python "${PACKAGE}/09_aggregate_genetic_context.py"
set +u
micromamba activate rnaseq
set -u
Rscript "${PACKAGE}/13_render_genetic.R"
set +u
micromamba activate snapatac2
set -u
python "${PACKAGE}/11_validate_genetics.py"
