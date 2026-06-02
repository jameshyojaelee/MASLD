#!/bin/bash
#SBATCH --job-name=prep_1kg_eur
#SBATCH --partition=cpu
#SBATCH --array=1-22
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/prep_1kg_eur_%a_%j.out
#SBATCH --error=RNA-seq/logs/prep_1kg_eur_%a_%j.err

# =============================================================================
# Convert 1000 Genomes Phase 1 VCFs to PLINK2 pgen format for EUR samples
# One chromosome per array task. Output: data/1kg_eur/chr{N}_eur.{pgen,psam,pvar}
# =============================================================================

set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID}
echo "=== Processing chromosome ${CHR} ==="
echo "Start: $(date)"

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
KG_DIR="/gpfs/commons/projects/1000Genomes_Phase1"
OUT_DIR="${BASE}/data/1kg_eur"
PANEL="${KG_DIR}/integrated_call_samples.20101123.ALL.panel"

mkdir -p "${OUT_DIR}" "${BASE}/RNA-seq/logs"

module load PLINK/2.0a5.13

# Step 1: Extract EUR sample IDs and create PLINK keep file
EUR_SAMPLES="${OUT_DIR}/eur_samples.txt"
EUR_KEEP="${OUT_DIR}/eur_samples_keep.txt"
if [[ ! -f "${EUR_SAMPLES}" ]]; then
  awk '$3 == "EUR" { print $1 }' "${PANEL}" > "${EUR_SAMPLES}"
  echo "  Extracted $(wc -l < "${EUR_SAMPLES}") EUR samples"
fi
# PLINK2 --keep needs two-column file (FID IID); VCF samples use ID as both
if [[ ! -f "${EUR_KEEP}" ]]; then
  awk '{ print $1, $1 }' "${EUR_SAMPLES}" > "${EUR_KEEP}"
fi

# Step 2: Find VCF for this chromosome
VCF="${KG_DIR}/ALL.chr${CHR}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz"
if [[ ! -f "${VCF}" ]]; then
  echo "ERROR: VCF not found: ${VCF}"
  exit 1
fi

# Step 3: Convert VCF → PLINK2 pgen with EUR samples only
# - --snps-only: exclude indels/SVs (COLOC uses SNPs)
# - --max-alleles 2: biallelic only
# - --maf 0.01: exclude rare variants (insufficient LD information)
# - --set-all-var-ids: format as chr:pos:ref:alt for unambiguous matching
# - --keep: EUR samples only
echo "  Converting VCF to pgen (EUR only, biallelic SNPs, MAF >= 0.01)..."

OUT_PREFIX="${OUT_DIR}/chr${CHR}_eur"

# Skip if already completed
if [[ -f "${OUT_PREFIX}.pgen" ]] && [[ -f "${OUT_PREFIX}.pvar" ]] && [[ -f "${OUT_PREFIX}.psam" ]]; then
  echo "  Output already exists, skipping: ${OUT_PREFIX}.pgen"
  echo "Done: $(date)"
  exit 0
fi

# Note: PLINK/2.0a5.13 module installs binary as 'plink' not 'plink2'
plink \
  --vcf "${VCF}" \
  --keep "${EUR_KEEP}" \
  --snps-only just-acgt \
  --max-alleles 2 \
  --maf 0.01 \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 20 \
  --make-pgen \
  --threads ${SLURM_CPUS_PER_TASK} \
  --out "${OUT_PREFIX}"

echo "  Variants: $(wc -l < "${OUT_PREFIX}.pvar")"
echo "  Samples:  $(wc -l < "${OUT_PREFIX}.psam")"
echo "Done: $(date)"
