#!/bin/bash
#SBATCH --job-name=N4_within_person
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --account=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# N4. Usage: sbatch --output=<out>/slurm-%j.log run_n4.sbatch <out_dir> <N1_dir|none>
set -euo pipefail
out=$1; n1=${2:-none}
repo=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
export OMP_NUM_THREADS=8 R_LIBS_USER=
test -s "$out/PRESPEC.md"; test -s "$out/pair_type_table.tsv"
test ! -e "$out/summary_by_group.tsv"
cp "$repo/scripts/analysis/hypotheses_2026q4/N4_within_person_agreement/n4_within_person_agreement.R" "$out/n4_within_person_agreement.as_run.R"
cp "$repo/scripts/analysis/hypotheses_2026q4/N4_within_person_agreement/n4_null_eligibility.R" "$out/n4_null_eligibility.as_run.R"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla \
  "$repo/scripts/analysis/hypotheses_2026q4/N4_within_person_agreement/test_n4_null_eligibility.R" "$repo"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla \
  "$repo/scripts/analysis/hypotheses_2026q4/N4_within_person_agreement/n4_within_person_agreement.R" "$out" "$n1"
date -u > "$out/COMPLETE.txt"
