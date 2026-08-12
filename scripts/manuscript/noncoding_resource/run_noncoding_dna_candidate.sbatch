#!/bin/bash
#SBATCH --job-name=ncdna
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncdna_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncdna_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
OUTPUT="${NONCODING_DNA_OUTPUT:-${PROJECT_ROOT}/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_${SLURM_JOB_ID}}"
IDENTITY="${PROJECT_ROOT}/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11/annotation/GENCODE_v49/run_19749483/gene_identity.tsv"

test ! -e "${OUTPUT}"

python3 scripts/manuscript/noncoding_resource/build_noncoding_dna_candidate.py \
  --output-dir "${OUTPUT}" \
  --credible-sets GWAS/finemapping/results/credible_sets.csv \
  --consequence-annotation GWAS/finemapping/results/credible_set_variant_consequences.csv \
  --trait-registry GWAS/finemapping/config/gwas_trait_tier.tsv \
  --gene-identity "${IDENTITY}" \
  --variant-liftover GWAS/finemapping/results/gwas_atac/variant_overlap_summary.csv \
  --abc-context GWAS/finemapping/results/gwas_atac/abc_variant_to_gene.csv \
  --atac-context GWAS/finemapping/results/gwas_atac/gwas_atac_variant_annotation.csv
