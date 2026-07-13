#!/usr/bin/env bash
#SBATCH --job-name=metafor
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/splicing_meta_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/logs/splicing_meta_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT=${ROOT}
micromamba run -n rnaseq Rscript "${ROOT}/GWAS/finemapping/src/107_splicing_meta_analysis.R"
test -s "${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/production/meta_analysis/random_effects_intron_meta.tsv"
