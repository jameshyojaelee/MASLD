#!/bin/bash
# CG0: ENSGID-keyed MAGMA gene analysis (SNP-wise mean, g1000_eur, 0 kb window)
# for the B4 trait set. PDFF uses the position-corrected Pazoki file.
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
cd $ROOT
OUT=$1
mkdir -p $OUT
MAGMA=tools/magma/magma
LUT=tools/magma/work/chr_pos_to_rsid.tsv
ANNOT=tools/magma/work/pops/magma_annot_ensgid.genes.annot
SS=GWAS/finemapping/data/sumstats
# id | file | N (gwas_registry.tsv N_tot; Ghouse from tools/magma _N_table)
TRAITS="ukbb_alt|$SS/UKBB_ALT_reformatted_hg19.tsv|343850
ukbb_ast|$SS/UKBB_AST_reformatted_hg19.tsv|343850
mvp_alt_eur|$SS/MVP_ALT_EUR_reformatted_hg19.tsv|406803
mvp_ast_eur|$SS/MVP_AST_EUR_reformatted_hg19.tsv|389192
pazoki_pdff|$SS/Pazoki_PDFF_reformatted_hg19.tsv|44867
ghouse_cirrhosis|$SS/Ghouse_Cirrhosis_reformatted_hg19.tsv|6000"
sha256sum $LUT $ANNOT $(echo "$TRAITS" | cut -d'|' -f2) GWAS/MR_Data/Yengo_2018_BMI_GWAS.txt.gz > $OUT/input_sha256.txt

# chr:pos -> rsID in one awk pass per trait (columns chromosome, position, pval)
echo "$TRAITS" | while IFS='|' read id f n; do
  awk -F'\t' -v N=$n 'NR==FNR{if(FNR>1)l[$1]=$2; next}
    FNR==1{for(i=1;i<=NF;i++)h[$i]=i; print "SNP\tP\tN"; next}
    {c=$h["chromosome"]; sub(/^chr/,"",c); k=c":"$h["position"]; p=$h["pval"];
     if((k in l) && p>0 && p<=1) print l[k]"\t"p"\t"N}' $LUT $f > $OUT/$id.magma_input.tsv &
done
wait
zcat GWAS/MR_Data/Yengo_2018_BMI_GWAS.txt.gz | awk 'NR==1{print "SNP\tP\tN"; next} {p=$9+0; if(p>0 && p<=1) print $3"\t"$9"\t"$10}' > $OUT/yengo_bmi.magma_input.tsv
wc -l $OUT/*.magma_input.tsv

for id in ukbb_alt ukbb_ast mvp_alt_eur mvp_ast_eur pazoki_pdff ghouse_cirrhosis yengo_bmi; do
  $MAGMA --bfile tools/magma/ref/g1000_eur --pval $OUT/$id.magma_input.tsv use=SNP,P ncol=N \
    --gene-annot $ANNOT --out $OUT/$id > $OUT/$id.stdout 2>&1 &
done
wait
grep -H -E "valid SNP p-values|genes" $OUT/*.log | head -40
