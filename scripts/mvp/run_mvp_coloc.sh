#!/bin/bash
#SBATCH --job-name=susie
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --array=0-593
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/logs/coloc_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/logs/coloc_%A_%a.err

# run_mvp_coloc.sh — SuSiE-COLOC for all 27 MVP strata x 22 autosomes.
# Array index 0..(27*22-1): stratum = idx/22 (row in registry-aligned manifest), chr = idx%22 + 1.
# (571.8/571.81 dropped 2026-07-02; manifest trimmed 32->27 rows 2026-07-04.)
# Canonical finemapping-tier COLOC: src/06_susie_coloc.R, registry-driven.
# LD_PANEL=polyfun -> EUR uses polyfun_eur; AFR/AMR/EAS fall to 1kg_<anc> via
# get_ld_base_dir(). eQTL = Broadaway EUR (hg19). MVP sumstats were lifted hg38->hg19
# upstream by format_mvp_for_coloc.R, so 06's inputs here are all hg19 (no liftOver at COLOC).

set -euo pipefail
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
FM="$BASE/GWAS/finemapping"
MAN="$BASE/GWAS/MR_Data/MVP/mvp_manifest.tsv"

module load PLINK/2.0a5.13
set +u   # micromamba activate.d (binutils) references unbound vars under set -u
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

N_CHR=22
SID=$(( SLURM_ARRAY_TASK_ID / N_CHR ))      # 0-based stratum index
CHR=$(( SLURM_ARRAY_TASK_ID % N_CHR + 1 ))  # 1..22
ROW=$(( SID + 2 ))                           # +1 header, +1 to 1-index
NAME=$(sed -n "${ROW}p" "$MAN" | cut -f1)
[ -z "$NAME" ] && { echo "No stratum at SID=$SID — done"; exit 0; }
# guard: skip cleanly if the harmonized sumstats file is absent (avoids a hard
# fread() error that would silently drop the stratum)
[ -f "$FM/data/sumstats/${NAME}_reformatted_hg19.tsv" ] || { echo "missing sumstats for $NAME — skip"; exit 0; }

echo "=== COLOC $NAME chr$CHR (task $SLURM_ARRAY_TASK_ID) ==="
cd "$FM"
export LD_PANEL=polyfun
GWAS_NAME="$NAME" CHR_FILTER="$CHR" Rscript src/06_susie_coloc.R
echo "=== done $NAME chr$CHR ==="
