#!/bin/bash
#SBATCH --job-name=cellranger
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --output=Analysis/Perturbation/scripts/saunders_reprocess/logs/08_anndata_%j.log

# Build Saunders Perturb-Multi mouse hepatocyte AnnData (~190K confident cells × ~20K genes)
# from 32 cellranger h5 GEX matrices + per-cell sgRNA assignment table.
# Output is GEARS-ready (sparse raw counts; obs.target_gene = mouse symbol; obs.human_ortholog).

set -o pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

eval "$(micromamba shell hook --shell bash)"
micromamba activate -p .envs/perturbation_gears

set -eu  # defer -u until after activate (MKL_INTERFACE_LAYER unbound during conda activate.d hooks)

python Analysis/Perturbation/scripts/saunders_reprocess/08_build_anndata.py \
  --per-cell  data/perturbation/datasets/saunders2025/processed/per_cell_sgrna.tsv.gz \
  --h5-glob   'data/perturbation/datasets/saunders2025/raw/h5/GSM8478327_*.h5' \
  --ortholog-map data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz \
  --out       data/perturbation/datasets/saunders2025/processed/saunders_perturb_hep_confident.h5ad \
  --confident-only
