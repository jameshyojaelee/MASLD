#!/bin/bash
#SBATCH --job-name=9135
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH -c 8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Methods/r1104_relabel_20261015/slurm_9135_%j.log
export LC_ALL=C.UTF-8 LINEAGES=hepatocyte
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export TMPDIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Methods/r1104_relabel_20261015/tmp; mkdir -p $TMPDIR
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Methods/r1104_relabel_20261015/relabel_check_sel.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/check/Methods/r1104_relabel_20261015/out_Cutf8_hep_backup 8
