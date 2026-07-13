#!/usr/bin/env bash
#SBATCH --job-name=saturation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/100_saturation_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/100_saturation_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export PYTHONNOUSERSITE=1 MASLD_PROJECT_ROOT="$ROOT"
mkdir -p "$ROOT/GWAS/finemapping/results/seqfunc/logs"
micromamba run -n rnaseq Rscript "$ROOT/GWAS/finemapping/src/99_haplotype_saturation_design.R"
micromamba run -n rnaseq python "$ROOT/GWAS/finemapping/src/100_saturation_sequence_prepare.py" \
  --manifest "$ROOT/GWAS/finemapping/results/seqfunc/haplotype_saturation/saturation_candidate_manifest.tsv" \
  --fasta /gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
  --out-dir "$ROOT/GWAS/finemapping/results/seqfunc/haplotype_saturation/sequences" \
  --n-targets "${N_TARGETS:-500}" --shard-size "${SHARD_SIZE:-25000}"
python3 "$ROOT/GWAS/finemapping/src/102a_repair_saturation_ids.py"
