#!/bin/bash
#SBATCH --job-name=Rscript
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --array=1-32
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/logs/harmonize_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/MVP/logs/harmonize_%A_%a.err

# run_mvp_harmonize.sh — array over mvp_manifest.tsv rows; reformat each MVP
# stratum into the finemapping hg19 COLOC schema. Requires mvp_manifest.tsv
# (produced by build_mvp_registry.R).

set -euo pipefail
BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
MAN="$BASE/GWAS/MR_Data/MVP/mvp_manifest.tsv"

set +u   # micromamba activate.d (binutils) references unbound vars under set -u
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

ROW=$((SLURM_ARRAY_TASK_ID + 1))           # +1 to skip header
LINE=$(sed -n "${ROW}p" "$MAN")
[ -z "$LINE" ] && { echo "No manifest row $ROW — done"; exit 0; }

NAME=$(echo "$LINE"   | cut -f1)
IN_GZ=$(echo "$LINE"  | cut -f2)
OUT_REL=$(echo "$LINE"| cut -f3)
OUT_ABS="$BASE/GWAS/finemapping/$OUT_REL"
mkdir -p "$(dirname "$OUT_ABS")"

echo "=== harmonize $NAME ==="
Rscript "$BASE/GWAS/MR_Data/MVP/format_mvp_for_coloc.R" "$IN_GZ" "$OUT_ABS"
echo "=== done $NAME -> $OUT_ABS ==="
