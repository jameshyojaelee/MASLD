#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# B-Ancestry: 1kGP reference PCs at the identity SNVs, projection of every MASLD library, and
# per-superpopulation tag REF/ALT heterozygote frequency and tag-lead LD. Runs anc0 -> anc3 in order.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j_<n>_ancestry_<mode>.out run_ancestry.sbatch smoke|full
#   smoke: chr20-22 -> ancestry/smoke-<UTC>/ (summary) and restricted ancestry/smoke-<UTC>/
#   full : chr1-22  -> ancestry/full-<UTC>/ (summary) and restricted ancestry/full-<UTC>/
#   (the first full run, job 21974425, wrote ancestry/ and restricted ancestry/work/; it is kept as is)
# Refuses to write into a summary directory that already holds reference_pca_summary.json.
set -euo pipefail
mode=$1
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
here=$repo/scripts/analysis/in_biopsy_allelic/ancestry
ex=$repo/Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-obs-20260926T155651Z
rs=$repo/RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted
a1=$repo/Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z
tagdir=$repo/Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-tags-20260925T194248Z
kgp=$repo/data/external/allelic_refs/kgp_highcov
ped=$kgp/20130606_g1k_3202_samples_ped_population.txt
py=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
threads=${SLURM_CPUS_PER_TASK:-8}
export PYTHONNOUSERSITE=1
module load bcftools/1.21 htslib/1.21 plink/2.0a5.13 >/dev/null 2>&1

ts=$(date -u +%Y%m%dT%H%M%SZ)
case "$mode" in
  smoke) outs=$ex/ancestry/smoke-$ts; outr=$rs/obs-20260926T155651Z/ancestry/smoke-$ts
         chroms=chr20,chr21,chr22 ;;
  full)  outs=$ex/ancestry/full-$ts; outr=$rs/obs-20260926T155651Z/ancestry/full-$ts
         chroms=$(seq -s, -f 'chr%g' 1 22) ;;
  *) echo "mode must be smoke or full" >&2; exit 2 ;;
esac
if [ -e "$outs/reference_pca_summary.json" ]; then echo "refusing to overwrite $outs" >&2; exit 3; fi
mkdir -p "$outs" "$outr"

"$py" "$here/anc0_regions.py" --targets "$a1/targets/targets.tsv.gz" --tags "$tagdir/t3/tags_final.tsv" \
  --ped "$ped" --chroms "$chroms" --out "$outr/sites"

bash "$here/anc1_extract_kgp.sh" "$outr/sites" "$kgp" "$chroms" "$outr/kgp_founders_sites.vcf.gz" "$threads" \
  "$outr/kgp_founders_tagspan.vcf.gz"

"$py" "$here/anc2_pca_project.py" --ref-vcf "$outr/kgp_founders_sites.vcf.gz" \
  --identity-sites "$outr/sites/identity_sites.tsv" --exclude-bed "$tagdir/excluded_spans_mhc_ig_tr.bed" \
  --tag-positions "$tagdir/t3/tags_final.tsv" "$tagdir/t3/tags_mappability.tsv" \
  --masld-vcf "$rs/a1-full-20260924T025858Z/vcf/identity_panel.merged.vcf.gz" \
  --crosswalk "$a1/identity/sample_to_individual.tsv" \
  --frozen-crosswalk "$a1/seal-20260929T042934Z/frozen_crosswalk.tsv" \
  --library-qc "$rs/obs-20260926T155651Z/qc/library_qc.tsv" --ped "$ped" \
  --work "$outr" --summary "$outs" --threads "$threads"

"$py" "$here/anc3_tag_ld_2pq.py" --ref-vcf "$outr/kgp_founders_sites.vcf.gz" \
  --tagspan-vcf "$outr/kgp_founders_tagspan.vcf.gz" --tags "$tagdir/t3/tags_final.tsv" \
  --ped "$ped" --unrelated "$outr/kgp_unrelated.keep" --exclude-bed "$tagdir/excluded_spans_mhc_ig_tr.bed" \
  --out "$outs"

{ module list 2>&1; bcftools --version | head -2; plink --version; "$py" -m pip freeze; } > "$outs/environment.txt"
echo "B-Ancestry $mode done: $outs $outr"
