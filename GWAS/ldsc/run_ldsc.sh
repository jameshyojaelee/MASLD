#!/bin/bash
#SBATCH --job-name=ldsc
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=GWAS/ldsc/logs/ldsc_%j.out
#SBATCH --error=GWAS/ldsc/logs/ldsc_%j.err
# Cross-trait LDSC genetic correlation (rg) + h2 across EUR liver-disease GWAS,
# for the GWAS genetic-architecture composite. Matrix = rg; lollipop uses COLOC
# counts (built separately). h2 is best-effort (FinnGen case counts approximate).
set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
mkdir -p GWAS/ldsc/logs GWAS/ldsc/munged GWAS/ldsc/rg
REF=GWAS/ldsc_ref/eur_w_ld_chr
MG=GWAS/ldsc/munged
run(){ micromamba run -n ldsc "$@"; }

# ── munge each trait: munge_one <name> <file> <snp> <a1(effect)> <a2> <p> <Nargs...> ──
munge_one(){
  local name=$1 file=$2 snp=$3 a1=$4 a2=$5 p=$6; shift 6
  [[ -f $MG/$name.sumstats.gz ]] && { echo "skip $name (munged)"; return; }
  local src=$file
  if zcat "$file" | head -1 | grep -q '^#chrom'; then      # FinnGen: strip leading #
    src=GWAS/ldsc/munged/${name}_src.tsv.gz
    zcat "$file" | sed '1s/^#//' | gzip > "$src"
  fi
  echo ">> munge $name"
  run munge_sumstats.py --sumstats "$src" --snp "$snp" --a1 "$a1" --a2 "$a2" \
      --signed-sumstats beta,0 --p "$p" "$@" --chunksize 500000 --out "$MG/$name" || echo "FAIL munge $name"
}

FG=GWAS/MR_Data/FinnGen; MR=GWAS/MR_Data
# NAFLD
munge_one ghodsian_nafld   $MR/Ghodsian_2021_NAFLD_harmonised.tsv.gz        hm_rsid hm_effect_allele hm_other_allele p_value --N-cas 8434 --N-con 770180
munge_one finngen_nafld    $FG/finngen_R12_NAFLD.gz                          rsids alt ref pval --N-cas 4900 --N-con 448000
# cirrhosis
munge_one finngen_cirr     $FG/finngen_R12_CHIRHEP_NAS.gz                    rsids alt ref pval --N-cas 2000 --N-con 450000
munge_one ghouse_cirr      $MR/Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz rsid effect_allele other_allele p_value --N-col samplesize
munge_one chen_cirr        $MR/Chen_2023_GOLDPlus_harmonised.tsv.gz          rsid effect_allele other_allele p_value --N 500000
# HCC
munge_one finngen_hcc      $FG/finngen_R12_C3_HEPATOCELLU_CARC_EXALLC.gz     rsids alt ref pval --N-cas 500 --N-con 452000
# liver enzymes (UKBB, continuous)
munge_one ukbb_alt         $MR/GCST90019492_UKBB_ALT_harmonised.tsv.gz       rsid effect_allele other_allele p_value --N 344000
munge_one ukbb_ast         $MR/GCST90019497_UKBB_AST_harmonised.tsv.gz       rsid effect_allele other_allele p_value --N 344000
munge_one ukbb_ggt         $MR/GCST90019507_UKBB_GGT_harmonised.tsv.gz       rsid effect_allele other_allele p_value --N 344000
# metabolic comorbidity
munge_one finngen_obesity  $FG/finngen_R12_E4_OBESITY.gz                     rsids alt ref pval --N-cas 40000 --N-con 410000

# ── per-trait h2 ──
TRAITS=(ghodsian_nafld finngen_nafld finngen_cirr ghouse_cirr chen_cirr finngen_hcc ukbb_alt ukbb_ast ukbb_ggt finngen_obesity)
AVAIL=(); for t in "${TRAITS[@]}"; do [[ -f $MG/$t.sumstats.gz ]] && AVAIL+=("$t"); done
echo "available munged: ${AVAIL[*]}"
n=${#AVAIL[@]}
for t in "${AVAIL[@]}"; do
  echo ">> h2 $t"
  run ldsc.py --h2 "$MG/$t.sumstats.gz" --ref-ld-chr $REF/ --w-ld-chr $REF/ \
      --out GWAS/ldsc/rg/h2_$t || echo "FAIL h2 $t"
done

# ── genetic correlation: anchor i vs all j>i (fills upper triangle) ──
for ((i=0; i<n-1; i++)); do
  anchor=${AVAIL[$i]}
  paths="$MG/${anchor}.sumstats.gz"
  for ((j=i+1; j<n; j++)); do paths+=",$MG/${AVAIL[$j]}.sumstats.gz"; done
  echo ">> rg anchored at $anchor"
  run ldsc.py --rg "$paths" --ref-ld-chr $REF/ --w-ld-chr $REF/ \
      --out GWAS/ldsc/rg/${anchor} || echo "FAIL rg $anchor"
done
echo "LDSC DONE"
