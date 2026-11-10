#!/bin/bash
#SBATCH --account=innovation --qos=innovation --job-name=code-review
#SBATCH --partition=cpu --prefer=v6 --cpus-per-task=4 --mem=8G --time=03:00:00
#SBATCH --array=1-493%50
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/logs/count/%a_%A.log
# SSX-B intronic counts, one staged included SteatoSITE library per task (PRESPEC_SSX_B.md; plan approval covers this array).
set -euo pipefail
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/scripts/ssx_b_samples.txt)
bash /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/scripts/ssx_b_count_one.sh $S /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/out_b_counts
