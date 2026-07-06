#!/bin/bash -l
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --partition=cpu,io
#SBATCH --job-name=mvploci
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/mvploci_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/logs/mvploci_%j.err
#
# Identify within-MVP multi-ancestry shared loci (09b). Streams GW-sig subsets
# from ~30GB of MVP sumstats via awk (cached under results/susiex_mvp/gwsig_cache).

set -eo pipefail
FM_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping"
cd "${FM_DIR}"
mkdir -p logs
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u
echo "[09b_mvp] start $(date)"
Rscript src/09b_identify_shared_loci_mvp.R
echo "[09b_mvp] done  $(date)"
