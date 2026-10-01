#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
# Tag REF/ALT het frequency and tag-lead LD per 1kGP superpopulation for the T1 v2 tag list: the tag
# steps (anc0, anc1, anc3) of B-Ancestry's run_ancestry.sbatch, unchanged scripts, public 1kGP only.
# Usage: sbatch --job-name=<n> --output=<out>/logs/%j.out run_anc_tags_v2.sbatch <out> <tags.tsv>
set -euo pipefail
out=$1; tags=$2
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
here=$repo/scripts/analysis/in_biopsy_allelic/ancestry
a1=$repo/Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z
tagdir=$repo/Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-tags-20260925T194248Z
kgp=$repo/data/external/allelic_refs/kgp_highcov
ped=$kgp/20130606_g1k_3202_samples_ped_population.txt
unrel=$repo/RNA-seq/Human/Patient_Cohorts/results/rna_genotypes_restricted/obs-20260926T155651Z/ancestry/full-20260929T142807Z/kgp_unrelated.keep
py=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
module load bcftools/1.21 htslib/1.21 plink/2.0a5.13 >/dev/null 2>&1
chroms=$(seq -s, -f 'chr%g' 1 22)
mkdir -p "$out/work"
"$py" "$here/anc0_regions.py" --targets "$a1/targets/targets.tsv.gz" --tags "$tags" --ped "$ped" --chroms "$chroms" --out "$out/work/sites"
bash "$here/anc1_extract_kgp.sh" "$out/work/sites" "$kgp" "$chroms" "$out/work/kgp_founders_sites.vcf.gz" 8 "$out/work/kgp_founders_tagspan.vcf.gz"
"$py" "$here/anc3_tag_ld_2pq.py" --ref-vcf "$out/work/kgp_founders_sites.vcf.gz" --tagspan-vcf "$out/work/kgp_founders_tagspan.vcf.gz" \
  --tags "$tags" --ped "$ped" --unrelated "$unrel" --exclude-bed "$tagdir/excluded_spans_mhc_ig_tr.bed" --out "$out"
{ module list 2>&1; bcftools --version | head -2; "$py" -m pip freeze; } > "$out/environment.txt"
