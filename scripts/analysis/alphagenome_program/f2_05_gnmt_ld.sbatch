#!/usr/bin/env bash
#SBATCH --job-name=agp_f2-haplotype-ld
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f2-haplotype-ld-%j.out
# Which of the four GNMT arms is a haplotype that EXISTS in Europeans, measured on this project's own
# 1000G EUR panel (379 samples, hg19). At r2 = 1 the single-variant arms are sequences nobody carries,
# so the additivity question is about the model's response, not about a segregating haplotype.
set -euo pipefail
module load PLINK/1.9b_6.21-x86_64
OUT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_program/f2-haplotype-20260914T230433Z
KG=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/1kg_eur
W="$OUT/raw/gnmt_ld"
mkdir -p "$W"
# hg19 positions from the panel this package built: rs9471976 42920549, rs11752813 42928017,
# rs2296805 42928758, rs2296804 42931261
awk -F'\t' '$4==42920549||$4==42928017||$4==42928758||$4==42931261{print $2}' "$KG/chr6_eur.bim" > "$W/gnmt_ids.txt"
cat "$W/gnmt_ids.txt"
plink --bfile "$KG/chr6_eur" --extract "$W/gnmt_ids.txt" --freq --out "$W/gnmt"
plink --bfile "$KG/chr6_eur" --extract "$W/gnmt_ids.txt" --r2 dprime --ld-window 10 --ld-window-kb 100 \
      --ld-window-r2 0 --out "$W/gnmt"
for A in $(cat "$W/gnmt_ids.txt"); do for B in $(cat "$W/gnmt_ids.txt"); do
  if [ "$A" != "$B" ]; then plink --bfile "$KG/chr6_eur" --ld "$A" "$B" --out "$W/ld_${A//:/_}__${B//:/_}" || true; fi
done; done
{ echo "== allele frequencies"; cat "$W/gnmt.frq"; echo; echo "== pairwise r2 / D'"; cat "$W/gnmt.ld"; echo;
  echo "== haplotype counts per pair"; grep -h -A12 "haplotype frequency" "$W"/ld_*.log || true; } > "$OUT/tables/gnmt_1kg_eur_ld.txt"
cat "$OUT/tables/gnmt_1kg_eur_ld.txt"
echo F2_LD_DONE
