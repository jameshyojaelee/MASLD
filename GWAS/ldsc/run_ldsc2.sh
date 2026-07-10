#!/bin/bash
#SBATCH --job-name=ldsc
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=GWAS/ldsc/logs/ldsc2_%j.out
#SBATCH --error=GWAS/ldsc/logs/ldsc2_%j.err
# Well-powered standard-EUR expansion: PDFF + MVP(ALT/AST/NAFLD/Cirrhosis, rsID-annotated
# from the HapMap3 chr:bp map). rg over the 9-trait set (FinnGen dropped — SNP-set mismatch).
set -euo pipefail
cd "${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
REF=GWAS/ldsc_ref/eur_w_ld_chr; MG=GWAS/ldsc/munged; MAP=GWAS/ldsc_ref/hm3_chrbp_to_rsid.tsv
SUM=GWAS/finemapping/data/sumstats
run(){ micromamba run -n ldsc "$@"; }

# ── annotate MVP (chr:pos hg19 -> rsID via HapMap3 map), keep mapped rows ──
annot_mvp(){  # <trait> <file>
  local key=$1 f=$2
  [[ -f $MG/mvp_${key}_src.tsv.gz ]] && return
  awk -F'\t' 'NR==FNR{m[$1]=$2;next} FNR==1{print "SNP\tA1\tA2\tBETA\tSE\tP";next}
    {k=$1":"$2; if(k in m) print m[k]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' "$MAP" "$f" | gzip > $MG/mvp_${key}_src.tsv.gz
  echo "mvp_$key annotated: $(zcat $MG/mvp_${key}_src.tsv.gz | wc -l) rows"
}
annot_mvp alt   $SUM/MVP_ALT_EUR_reformatted_hg19.tsv
annot_mvp ast   $SUM/MVP_AST_EUR_reformatted_hg19.tsv
annot_mvp nafld $SUM/MVP_NAFLD_EUR_reformatted_hg19.tsv
annot_mvp cirr  $SUM/MVP_Cirrhosis_EUR_reformatted_hg19.tsv

munge_new(){  # <name> <src> <extra munge args...>
  local name=$1 src=$2; shift 2
  [[ -f $MG/$name.sumstats.gz ]] && { echo "skip $name"; return; }
  run munge_sumstats.py --sumstats "$src" --snp SNP --a1 A1 --a2 A2 --signed-sumstats BETA,0 --p P \
      "$@" --chunksize 500000 --out $MG/$name 2>&1 | grep -iE "remain|chi|Writing|error" | tail -2 || echo "FAIL $name"
}
munge_new mvp_alt   $MG/mvp_alt_src.tsv.gz   --N 200000
munge_new mvp_ast   $MG/mvp_ast_src.tsv.gz   --N 200000
munge_new mvp_nafld $MG/mvp_nafld_src.tsv.gz --N-cas 30000 --N-con 350000
munge_new mvp_cirr  $MG/mvp_cirr_src.tsv.gz  --N-cas 10000 --N-con 370000
# PDFF (rsID variant_id, per-SNP N)
[[ -f $MG/pdff.sumstats.gz ]] || run munge_sumstats.py --sumstats GWAS/MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz \
   --snp variant_id --a1 effect_allele --a2 other_allele --signed-sumstats beta,0 --p p_value --N-col n \
   --chunksize 500000 --out $MG/pdff 2>&1 | grep -iE "remain|chi|Writing|error" | tail -2

# ── h2 + rg over the well-powered 9-trait set ──
T=(ghodsian_nafld mvp_nafld pdff ghouse_cirr mvp_cirr ukbb_alt mvp_alt ukbb_ast mvp_ast)
A=(); for t in "${T[@]}"; do [[ -f $MG/$t.sumstats.gz ]] && A+=("$t"); done
echo "well-powered set: ${A[*]}"; n=${#A[@]}
for t in "${A[@]}"; do run ldsc.py --h2 $MG/$t.sumstats.gz --ref-ld-chr $REF/ --w-ld-chr $REF/ --out GWAS/ldsc/rg/h2_$t >/dev/null 2>&1 || echo "h2 FAIL $t"; done
for ((i=0;i<n-1;i++)); do p="$MG/${A[$i]}.sumstats.gz"; for ((j=i+1;j<n;j++)); do p+=",$MG/${A[$j]}.sumstats.gz"; done
  run ldsc.py --rg "$p" --ref-ld-chr $REF/ --w-ld-chr $REF/ --out GWAS/ldsc/rg/wp_${A[$i]} >/dev/null 2>&1 && echo "rg $i ok" || echo "rg $i FAIL"; done
echo "LDSC2 DONE"
