#!/bin/bash
# Ancestry step 1: 1kGP high-coverage phased hardcalls for the founders at the anc0 positions.
#
# Inputs : <sites_dir> from anc0 (founders.txt, regions.<chrom>.tsv, tag_regions.<chrom>.tsv), the
#          per-chromosome 1kGP high-coverage phased panel VCFs, a comma list of chromosomes.
# Outputs: <out_vcf> (bgzipped, indexed): founders only, biallelic SNV records whose POS is one
#          of the anc0 positions (--regions-overlap pos, so deletions spanning a site are not pulled
#          in). GT is left exactly as in the panel (phased hardcalls); nothing is recomputed.
#          <out_tagspan_vcf> (bgzipped, indexed): founders only, EVERY panel record (SNV, indel, SV,
#          split multiallelic) whose REF span covers a tag SNV position (--regions-overlap variant).
#          anc3 uses it to count haplotypes that carry another base or a deletion at the tag base.
# Rules  : reads no MASLD data. Needs bcftools on PATH.
# Usage  : anc1_extract_kgp.sh <sites_dir> <panel_dir> <chrom_list> <out_vcf> <threads> <out_tagspan_vcf>
set -euo pipefail
sites=$1; panel=$2; chroms=$3; out=$4; threads=$5; out_span=$6
tmp=$(dirname "$out")/kgp_by_chrom
mkdir -p "$tmp"

extract() {
  local c=$1
  local vcf="$panel/1kGP_high_coverage_Illumina.$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  bcftools view -R "$sites/regions.$c.tsv" --regions-overlap pos -S "$sites/founders.txt" \
    -v snps -m2 -M2 -Ob -o "$tmp/$c.bcf" "$vcf"
  bcftools index -f "$tmp/$c.bcf"
  if [ -s "$sites/tag_regions.$c.tsv" ]; then
    bcftools view -R "$sites/tag_regions.$c.tsv" --regions-overlap variant -S "$sites/founders.txt" \
      -Ob -o "$tmp/$c.tagspan.bcf" "$vcf"
  else  # no tag on this chromosome: header-only file keeps the concat simple
    bcftools view -h -S "$sites/founders.txt" -Ob -o "$tmp/$c.tagspan.bcf" "$vcf"
  fi
  bcftools index -f "$tmp/$c.tagspan.bcf"
}
export -f extract
export sites panel tmp

echo "$chroms" | tr ',' '\n' | xargs -P "$threads" -I{} bash -c 'extract {}'
bcftools concat --threads "$threads" -Oz -o "$out" $(echo "$chroms" | tr ',' '\n' | sed "s|.*|$tmp/&.bcf|")
bcftools index -t -f "$out"
bcftools concat --threads "$threads" -Oz -o "$out_span" $(echo "$chroms" | tr ',' '\n' | sed "s|.*|$tmp/&.tagspan.bcf|")
bcftools index -t -f "$out_span"
echo "anc1: $(bcftools view -H "$out" | wc -l) records, $(bcftools query -l "$out" | wc -l) samples -> $out"
echo "anc1: $(bcftools view -H "$out_span" | wc -l) tag-span records -> $out_span"
