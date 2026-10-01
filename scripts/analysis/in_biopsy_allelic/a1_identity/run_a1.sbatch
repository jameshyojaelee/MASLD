#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
# Job name and log path are set on the sbatch command line (numbers only).
#
# A1 identity panel for the in-biopsy allelic model. Reads no outcome.
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out run_a1.sbatch <smoke|full> <exec_dir> <restricted_dir>
#   smoke: 20 genes, 2 BAMs from each of GSE130970 and GSE193066
#   full : 6,000 genes, all 1,466 BAMs in 10 cohorts
set -euo pipefail

mode=$1
exec_dir=$2
restricted=$3
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
src=$repo/scripts/analysis/in_biopsy_allelic/a1_identity
py=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
export PLINK2=plink
module load bcftools/1.21 htslib/1.21 plink/2.0a5.13 >/dev/null 2>&1

if [ "$mode" = smoke ]; then
  cohorts=(GSE130970 GSE193066); n_genes=20; max_bams=2; min_sites=10; min_shared=10
else
  cohorts=(GSE126848 GSE130970 GSE135251 GSE162694 GSE167523 GSE174478 GSE193066 GSE213621 GSE240729 PRJNA512027)
  n_genes=6000; max_bams=0; min_sites=2000; min_shared=1000
fi

mkdir -p "$exec_dir" "$restricted/vcf" "$restricted/bai"
chmod 750 "$restricted"
{ module list 2>&1; "$py" -m pip freeze; bcftools --version | awk 'NR==1'; plink --version; } > "$exec_dir/environment.txt"

"$py" "$src/01_targets.py" \
  --counts "$repo/RNA-seq/Human/Patient_Cohorts/results/GSE130970/counts/featurecounts/gene_counts.txt" \
  --gtf "$repo/Analysis/MASLD_Model_Benchmark/executions/reference-build-21062075/gencode.v49.primary_analysis.annotation.gtf.gz" \
  --dbsnp /gpfs/commons/home/jameslee/reference_genome/dbsnp/GCF_000001405.40.gz \
  --out "$exec_dir/targets" --n-genes "$n_genes"

# Index every BAM once, 16 at a time, into the restricted index dir.
root=$repo/RNA-seq/Human/Patient_Cohorts/results
for c in "${cohorts[@]}"; do
  ls "$root/$c"/alignments/star/*/*.Aligned.sortedByCoord.out.bam | awk -v n="$max_bams" 'n == 0 || NR <= n' \
    | awk -v c="$c" -v d="$restricted/bai" '{n=split($0,p,"/"); r=p[n]; sub(/\.Aligned\.sortedByCoord\.out\.bam$/,"",r); print $0"\t"d"/"c"/"r".bam"}'
done > "$exec_dir/index_jobs.tsv"
cut -f2 "$exec_dir/index_jobs.tsv" | xargs -n1 dirname | sort -u | xargs mkdir -p
awk -F'\t' '{print $1"\n"$2}' "$exec_dir/index_jobs.tsv" \
  | xargs -n2 -P 16 sh -c 'ln -sfn "$1" "$2"; [ -s "$2.bai" ] || /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools index -o "$2.bai" "$2"' sh

printf '%s\n' "${cohorts[@]}" \
  | xargs -P 8 -I{} bash "$src/02_call_cohort.sh" {} "$exec_dir/targets" "$restricted/vcf" "$restricted/bai" "$max_bams"

"$py" "$src/03_identity.py" --vcf-dir "$restricted/vcf" --cohorts "${cohorts[@]}" \
  --gtf "$repo/Analysis/MASLD_Model_Benchmark/executions/reference-build-21062075/gencode.v49.primary_analysis.annotation.gtf.gz" \
  --metadata "$repo/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv" \
  --summary-dir "$exec_dir/identity" --min-sites "$min_sites" --min-shared-snps "$min_shared" --threads 16
echo "A1 $mode done"
