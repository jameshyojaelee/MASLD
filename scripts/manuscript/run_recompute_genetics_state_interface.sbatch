#!/bin/bash
#SBATCH --job-name=gen_state_treat
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/gen_state_treat_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/gen_state_treat_%j.err

set -euo pipefail
umask 027
: "${GENETICS_STATE_INTERFACE_OUT:?Set a fresh output directory}"
test ! -e "${GENETICS_STATE_INTERFACE_OUT}"
project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
"${rnaseq_env}/bin/python3" \
  "${project_root}/scripts/manuscript/recompute_genetics_state_interface.py" \
  --out "${GENETICS_STATE_INTERFACE_OUT}"
test -s "${GENETICS_STATE_INTERFACE_OUT}/matched_treat_by_arm.tsv"
test -s "${GENETICS_STATE_INTERFACE_OUT}/matched_pairs_by_arm.tsv"
test -s "${GENETICS_STATE_INTERFACE_OUT}/threshold_free_by_arm.tsv"
echo "GENETICS_STATE_TREAT_RECOMPUTE_COMPLETE"
