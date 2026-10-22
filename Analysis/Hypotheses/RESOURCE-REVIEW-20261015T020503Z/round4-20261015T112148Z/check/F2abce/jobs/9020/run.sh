#!/bin/bash
set -euo pipefail
CK=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/F2abce
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 TMPDIR=$CK/tmp
OUT=$CK/out
echo "start $(date -Is) on $(hostname)"
$ENV/python $CK/scripts/01_prep.py $CK
# stream-filter the two fine-mapping tables to keys present in main COLOC rows
FM=$ROOT/GWAS/finemapping/results/combined_finemapping.csv
FT=$ROOT/GWAS/finemapping/results/combined_finemapping_coloc_targeted.csv
head -1 $FM > $OUT/fm_canon_header.txt
head -1 $FT > $OUT/fm_targ_header.txt
# canonical: chromosome 1, position 2, study 7 ; NF must be 27
mawk -F, 'NR==FNR{k[$1","$2","$3]=1;next} FNR==1{print;next} NF!=27{bad++; next} (($7","$1","$2) in k){print} END{print "bad_nf_canon",bad+0 > "/dev/stderr"}' $OUT/keys.csv $FM > $OUT/fm_canon_filtered.csv
# targeted: chromosome 1, position 2, study 6 ; NF must be 21
mawk -F, 'NR==FNR{k[$1","$2","$3]=1;next} FNR==1{print;next} NF!=21{bad++; next} (($6","$1","$2) in k){print} END{print "bad_nf_targ",bad+0 > "/dev/stderr"}' $OUT/keys.csv $FT > $OUT/fm_targ_filtered.csv
wc -l $OUT/fm_canon_filtered.csv $OUT/fm_targ_filtered.csv
# R16-08: r2 rows for the timed-out stratum, chr1
R2=$ROOT/GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv
mawk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $1=="2023_36280732_NAFLD_UKBB_EUR"{n[$4]++; if($(h["PP.H4.susie"])!="" && $(h["PP.H4.susie"])!="NA"){s[$4]++; if($(h["PP.H4.susie"])+0>0.5) p[$4]++}} $1=="2023_36280732_NAFLD_Intermountain_EUR"&&$4=="1"{im++} END{print "UKBB_EUR chr1 rows",n["1"],"susie_nonmissing",s["1"]+0,"pp4susie_gt_0.5",p["1"]+0; print "Intermountain chr1 rows", im}' $R2 > $OUT/r16_08_r2_counts.txt
md5sum $ROOT/GWAS/finemapping/results/susie_coloc_rerun/2023_36280732_NAFLD_UKBB_EUR/susie_coloc_chr1.csv $ROOT/GWAS/finemapping/results/susie_coloc/2023_36280732_NAFLD_UKBB_EUR/susie_coloc_chr1.csv > $OUT/r16_08_md5.txt
echo "done $(date -Is)"
