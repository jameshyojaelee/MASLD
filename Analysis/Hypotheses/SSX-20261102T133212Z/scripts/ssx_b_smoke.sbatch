#!/bin/bash
#SBATCH --account=innovation --qos=innovation --job-name=code-review
#SBATCH --partition=cpu --prefer=v6 --cpus-per-task=4 --mem=8G --time=03:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/logs/ssx_b_smoke_%j.log
set -euo pipefail
for S in S0001 $(sed -n 200p /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/scripts/ssx_b_samples.txt); do
  bash /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/scripts/ssx_b_count_one.sh $S /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/out_b_smoke; done
echo EXIT_OK
