#!/usr/bin/env bash
#SBATCH --job-name=leafcutter
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/annotate_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/logs/annotate_%j.err
set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export MASLD_PROJECT_ROOT=${ROOT}
micromamba run -n rnaseq Rscript "${ROOT}/GWAS/finemapping/src/seqfunc_v2/splicing/07_annotate_replicate.R"
test -s "${ROOT}/GWAS/finemapping/results/seqfunc/disease_splicing/joint_v2/primary_events_annotated.tsv.gz"
