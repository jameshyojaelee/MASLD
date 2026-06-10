#!/bin/bash
#SBATCH --job-name=figrender
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:40:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/figrender_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/figures/logs/figrender_%j.err

set -eo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p "$BASE/scripts/figures/logs"
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "######## figS_ccc_shape_clusters.R ########"
Rscript "$BASE/scripts/figures/figS_ccc_shape_clusters.R"

echo "######## fig_rora_case_study.R ########"
Rscript "$BASE/scripts/figures/fig_rora_case_study.R"
