#!/bin/bash
#SBATCH --job-name=ldsc
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=GWAS/ldsc/logs/ldsc3_%j.out
#SBATCH --error=GWAS/ldsc/logs/ldsc3_%j.err
# Re-munge ALL traits WITH --merge-alleles (HapMap3 allele reference) to fix the
# incompatible-allele failures (MVP rg NA + PDFF sign). Then h2 + rg over the full set.
set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
REF=GWAS/ldsc_ref/eur_w_ld_chr; MG=GWAS/ldsc/munged; MAP=GWAS/ldsc_ref/hm3_chrbp_to_rsid.tsv
FG=GWAS/MR_Data/FinnGen; MR=GWAS/MR_Data; SUM=GWAS/finemapping/data/sumstats
MA=GWAS/ldsc_ref/hm3_merge_alleles.txt
run(){ micromamba run -n ldsc "$@"; }

# ── build HapMap3 merge-alleles (SNP A1 A2) from UKBB ALT alleles ∩ HapMap3 ──
if [[ ! -f $MA ]]; then
  echo "SNP	A1	A2" > $MA
  zcat $MR/GCST90019492_UKBB_ALT_harmonised.tsv.gz | \
    awk -F'\t' 'NR==FNR{hm[$2]=1;next} FNR==1{for(i=1;i<=NF;i++)c[$i]=i;next} ($c["rsid"] in hm){print $c["rsid"]"\t"$c["effect_allele"]"\t"$c["other_allele"]}' \
    $MAP - >> $MA
  echo "merge-alleles rows: $(wc -l < $MA)"
fi

rm -f $MG/*.sumstats.gz     # force re-munge with --merge-alleles

munge(){  # <name> <file> <snp> <a1> <a2> <p> <Nargs...>
  local name=$1 file=$2 snp=$3 a1=$4 a2=$5 p=$6; shift 6
  local src=$file
  if zcat "$file" 2>/dev/null | head -1 | grep -q '^#chrom'; then
    src=$MG/${name}_src.tsv.gz; zcat "$file" | sed '1s/^#//' | gzip > "$src"; fi
  run munge_sumstats.py --sumstats "$src" --snp "$snp" --a1 "$a1" --a2 "$a2" --signed-sumstats beta,0 \
      --p "$p" "$@" --merge-alleles $MA --chunksize 500000 --out $MG/$name 2>&1 | grep -iE "remain|Writing|error" | tail -1 || echo "FAIL $name"
}
mungesrc(){  # pre-cleaned SNP/A1/A2/BETA/P src
  local name=$1 src=$2; shift 2
  run munge_sumstats.py --sumstats "$src" --snp SNP --a1 A1 --a2 A2 --signed-sumstats BETA,0 --p P \
      "$@" --merge-alleles $MA --chunksize 500000 --out $MG/$name 2>&1 | grep -iE "remain|Writing|error" | tail -1 || echo "FAIL $name"
}

mungesrc ghodsian_nafld $MG/ghodsian_src.tsv.gz --N-cas 8434 --N-con 770180
mungesrc mvp_alt   $MG/mvp_alt_src.tsv.gz   --N 200000
mungesrc mvp_ast   $MG/mvp_ast_src.tsv.gz   --N 200000
mungesrc mvp_nafld $MG/mvp_nafld_src.tsv.gz --N-cas 30000 --N-con 350000
mungesrc mvp_cirr  $MG/mvp_cirr_src.tsv.gz  --N-cas 10000 --N-con 370000
munge ukbb_alt     $MR/GCST90019492_UKBB_ALT_harmonised.tsv.gz rsid effect_allele other_allele p_value --N 344000
munge ukbb_ast     $MR/GCST90019497_UKBB_AST_harmonised.tsv.gz rsid effect_allele other_allele p_value --N 344000
munge ghouse_cirr  $MR/Ghouse_Cirrhosis/GCST90319877_harmonised_hg38.tsv.gz rsid effect_allele other_allele p_value --N-col samplesize
munge pdff         $MR/GCST90267352_PDFF_Pazoki2022.tsv.gz variant_id effect_allele other_allele p_value --N-col n
munge finngen_nafld $FG/finngen_R12_NAFLD.gz rsids alt ref pval --N-cas 4900 --N-con 448000
munge finngen_cirr  $FG/finngen_R12_CHIRHEP_NAS.gz rsids alt ref pval --N-cas 2000 --N-con 450000
munge finngen_obesity $FG/finngen_R12_E4_OBESITY.gz rsids alt ref pval --N-cas 40000 --N-con 410000

# ── h2 + rg over the full set ──
T=(ghodsian_nafld mvp_nafld finngen_nafld pdff ghouse_cirr mvp_cirr finngen_cirr ukbb_alt mvp_alt ukbb_ast mvp_ast finngen_obesity)
A=(); for t in "${T[@]}"; do [[ -f $MG/$t.sumstats.gz ]] && A+=("$t"); done
echo "final set: ${A[*]}"; n=${#A[@]}
for t in "${A[@]}"; do run ldsc.py --h2 $MG/$t.sumstats.gz --ref-ld-chr $REF/ --w-ld-chr $REF/ --out GWAS/ldsc/rg/h2_$t >/dev/null 2>&1 || echo "h2 FAIL $t"; done
for ((i=0;i<n-1;i++)); do p="$MG/${A[$i]}.sumstats.gz"; for ((j=i+1;j<n;j++)); do p+=",$MG/${A[$j]}.sumstats.gz"; done
  run ldsc.py --rg "$p" --ref-ld-chr $REF/ --w-ld-chr $REF/ --out GWAS/ldsc/rg/ma_${A[$i]} >/dev/null 2>&1 && echo "rg $i ok" || echo "rg $i FAIL"; done
echo "LDSC3 DONE"
