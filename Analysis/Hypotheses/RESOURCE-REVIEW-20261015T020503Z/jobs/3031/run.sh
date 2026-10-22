#!/bin/bash
# gate-L06-02: stream GWAS summary files and keep rows at 3,025 chr:pos keys
# (UKBB ALT/AST 500 kb clumped leads, PNPLA3 rs738409, TM6SF2 rs58542926, seeded null SNVs).
set -euo pipefail
J=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/jobs/3031
G=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/data
date
filt() { # $1 = source path, $2 = output name
  mawk -F'\t' 'NR==FNR{k[$1":"$2]=1; next} FNR==1{print; next} ($1":"$2) in k {print}' "$J/keys.tsv" "$1" > "$J/out/$2.tsv"
  echo "$2 $(wc -l < "$J/out/$2.tsv")"
}
export -f filt; export J
{
for t in ALT AST; do
  for a in EUR AFR AMR EAS; do echo "$G/sumstats_af/MVP_${t}_${a}_af.tsv MVP_${t}_${a}"; done
  echo "$G/sumstats_af/BBJ_${t}_af.tsv BBJ_${t}"
  echo "$G/sumstats_af/PanUKBB_AFR_${t}_af.tsv PanUKBB_AFR_${t}"
  echo "$G/sumstats_af/PanUKBB_CSA_${t}_af.tsv PanUKBB_CSA_${t}"
  echo "$G/sumstats/UKBB_${t}_reformatted_hg19.tsv UKBB_${t}"
done
echo "$G/sumstats/MVP_ALT_EUR_reformatted_hg19.tsv MVP_ALT_EUR_reformatted"
} | xargs -P 4 -L 1 bash -c 'filt "$0" "$1"'
date
