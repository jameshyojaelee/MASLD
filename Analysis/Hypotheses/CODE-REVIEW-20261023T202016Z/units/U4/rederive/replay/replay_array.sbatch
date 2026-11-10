#!/bin/bash
#SBATCH --job-name=u4_replay
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/logs/%A_%a.out
# --mem and --array are set on the command line by submit_array.sh from the
# single-task timing of job 22383429.
set -eo pipefail
/usr/bin/time -v /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/replay/replay_task.sh ${SLURM_ARRAY_TASK_ID}
