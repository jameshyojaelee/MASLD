#!/bin/bash
#SBATCH --job-name=rv7
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV7-20261102T145824Z/logs/%x_%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV7-20261102T145824Z
OUT=$RUN/out-$(date -u +%Y%m%dT%H%M%SZ)
mkdir "$OUT"
cd "$RUN" && sha256sum -c PRESPEC_RV7.md.sha256
sha256sum scripts/rv7_gnmt_conditional.R > "$OUT/script_sha256.txt"
module load plink/1.9
plink --bfile /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/data/ld_ref/1kg_eur/chr6_eur \
  --chr 6 --from-bp 41928548 --to-bp 43928108 --keep-allele-order --recode A --out "$OUT/ld_window"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/rv7_gnmt_conditional.R "$OUT"
