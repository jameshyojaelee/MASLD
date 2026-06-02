#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=cpu
#SBATCH --job-name=agg4way
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/agg4way_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/agg4way_%j.err
#
# Aggregate SuSiEX-4way + MESuSiE-4way into per-gene/per-variant CSVs and
# render the 2-way vs 4-way comparison figure.
#
# Usage:
#   sbatch --dependency=afterany:<sx_array>:<ms_array> src/run_aggregate_4way.sh

set -eo pipefail
FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
PROJ="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${PROJ}"

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

echo "[agg4way] start $(date)"
echo "[agg4way] aggregating SuSiEX 4-way..."
Rscript GWAS/finemapping/src/11_aggregate_susiex_4way.R
echo "[agg4way] aggregating MESuSiE 4-way..."
Rscript GWAS/finemapping/src/11b_aggregate_mesusie_4way.R
echo "[agg4way] rendering 2-way vs 4-way figure..."
Rscript scripts/figures/presentation_4way_ancestry.R
echo "[agg4way] done $(date)"
