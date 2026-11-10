#!/bin/bash
#SBATCH --job-name=gen-LD-pilot
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/genetics/pilot-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/genetics/pilot-%j.err
set -euo pipefail
work=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/genetics
cd "$work"
sha256sum -c SPEC.sha256
sha256sum -c pilot_sources.sha256
sha256sum -c pilot_code.sha256
export PATH=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin:$PATH
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export EQTL_SUSIE_DIR=$MASLD_PROJECT_ROOT/GWAS/finemapping/results/eqtl_susie_polyfun
export LD_PANEL=polyfun
unset UKBB_LD_DIR
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
export COLOC_GENE_FILTER=$work/pilot_targets.tsv
pilot_output=$work/pilot_${SLURM_JOB_ID}
mkdir "$pilot_output"
export COLOC_OUT_ROOT=$pilot_output/coloc
export GENETICS_PILOT_EXPORT_DIR=$pilot_output/exports
Rscript -e 'library(coloc); library(susieR); library(data.table); sessionInfo()' > "$pilot_output/sessionInfo.txt"
tail -n +2 "$work/pilot_targets.tsv" | while IFS=$'\t' read -r study chromosome ensembl gene role; do
  echo "Reconstructing $study $gene chr$chromosome"
  GWAS_NAME="$study" CHR_FILTER="$chromosome" Rscript "$work/pilot_runner_snapshot.R"
done
python "$work/check_pilot.py" "$pilot_output"
sha256sum -c pilot_sources.sha256
