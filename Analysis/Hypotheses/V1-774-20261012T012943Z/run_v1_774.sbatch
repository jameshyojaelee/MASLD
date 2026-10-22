#!/bin/bash
#SBATCH --job-name=v1_774_gc
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/V1-774-20261012T012943Z/slurm-%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export V1_OUT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V1-774-20261012T012943Z
export V1_RES=$V1_OUT/results/run-$(date -u +%Y%m%dT%H%M%SZ)
mkdir "$V1_RES"
cp "$V1_OUT/v1_774_gc_covariate.R" "$V1_RES/v1_774_gc_covariate.as_run.R"
sha256sum "$V1_OUT/v1_774_gc_covariate.R" "$V1_OUT/PRESPEC.md" "$V1_OUT/inputs/catalog_r8_concordant_genes.tsv" > "$V1_RES/script_sha256.txt"
echo "job $SLURM_JOB_ID host $(hostname) start $(date -u)"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$V1_OUT/v1_774_gc_covariate.R"
echo "end $(date -u)"
