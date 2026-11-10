#!/bin/bash
#SBATCH --job-name=u4_neg_recheck
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --array=1-345%10
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/neg_recheck_20261024T111335Z/logs/%A_%a.out
set -eo pipefail
sha256sum /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/src/06_susie_coloc.R /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/GWAS/finemapping/src/coloc_signal_eligibility.R
/usr/bin/time -v bash /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/neg_recheck_20261024T111335Z/task.sh ${SLURM_ARRAY_TASK_ID}
