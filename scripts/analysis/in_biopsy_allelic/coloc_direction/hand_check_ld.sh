#!/bin/bash
# Independent re-derivation of one LD r from the 1kGP phased panel (bcftools + awk, no Python code shared
# with coloc_direction.py). Args: chrom posA alleleA posB alleleB superpop. Prints r between alleleA at posA
# and alleleB at posB over founders of the superpopulation (both haplotypes of each founder).
set -euo pipefail
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
K=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/allelic_refs/kgp_highcov
chrom=$1; pa=$2; aa=$3; pb=$4; ab=$5; sp=$6
samples=$(awk -v sp="$sp" 'NR>1 && $7==sp && $3=="0" && $4=="0" {print $2}' $K/20130606_g1k_3202_samples_ped_population.txt | paste -sd,)
vcf=$K/1kGP_high_coverage_Illumina.$chrom.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
bcftools query -r "$chrom:$pa,$chrom:$pb" -s "$samples" -i 'N_ALT=1' -f '%POS\t%REF\t%ALT[\t%GT]\n' "$vcf" |
awk -v pa=$pa -v pb=$pb -v aa=$aa -v ab=$ab 'BEGIN{FS="\t"}
  $1==pa && ($2==aa || $3==aa) {na++; for(i=4;i<=NF;i++){split($i,g,"|"); x[2*i]=(g[1]==1)?($3==aa):($2==aa); x[2*i+1]=(g[2]==1)?($3==aa):($2==aa)}; print "A", $1, $2, $3}
  $1==pb && ($2==ab || $3==ab) {nb++; for(i=4;i<=NF;i++){split($i,g,"|"); y[2*i]=(g[1]==1)?($3==ab):($2==ab); y[2*i+1]=(g[2]==1)?($3==ab):($2==ab)}; print "B", $1, $2, $3}
  END{ if(na!=1||nb!=1){print "records", na, nb; exit 1}
       for(k in x){n++; sx+=x[k]; sy+=y[k]; sxy+=x[k]*y[k]; sxx+=x[k]*x[k]; syy+=y[k]*y[k]}
       r=(n*sxy-sx*sy)/sqrt((n*sxx-sx*sx)*(n*syy-sy*sy)); printf "haplotypes %d freqA %.4f freqB %.4f r %.6f\n", n, sx/n, sy/n, r }'
