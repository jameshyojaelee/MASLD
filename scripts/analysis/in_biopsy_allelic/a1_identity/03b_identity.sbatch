#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
# A1 step 3 after the chunk array: concatenate chunks per cohort, then 03_identity.py.
# Usage: sbatch --job-name=<n> --dependency=afterok:<array> --output=<exec>/logs/%j.out 03b_identity.sbatch <exec_dir> <restricted_dir>
set -euo pipefail
ex=$1; rs=$2
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
module load bcftools/1.21 htslib/1.21 plink/2.0a5.13 >/dev/null 2>&1
export PYTHONNOUSERSITE=1 PLINK2=plink
cohorts=(GSE126848 GSE130970 GSE135251 GSE162694 GSE167523 GSE174478 GSE193066 GSE213621 GSE240729 PRJNA512027)
mkdir -p "$rs/vcf"
for c in "${cohorts[@]}"; do
  bcftools concat -Oz -o "$rs/vcf/$c.vcf.gz" $(for k in $(seq 0 9); do printf '%s ' "$rs/vcf_chunks/$c.chunk$k.vcf.gz"; done)
  bcftools index -t "$rs/vcf/$c.vcf.gz"
done
{ module list 2>&1; /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze; } > "$ex/environment.txt"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python "$repo/scripts/analysis/in_biopsy_allelic/a1_identity/03_identity.py" \
  --vcf-dir "$rs/vcf" --cohorts "${cohorts[@]}" \
  --gtf "$repo/Analysis/MASLD_Model_Benchmark/executions/reference-build-21062075/gencode.v49.primary_analysis.annotation.gtf.gz" \
  --metadata "$repo/RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv" \
  --summary-dir "$ex/identity" --min-sites 1000 --min-shared-snps 500 --threads 8
echo "A1 identity done"
