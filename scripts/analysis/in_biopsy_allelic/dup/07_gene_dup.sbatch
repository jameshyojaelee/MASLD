#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
# Job name, array range and log path are set on the sbatch command line (numbers only).
#
# B-Dup step 7, one library per array task (task = the "task" column of libraries.tsv).
# Runs 07_gene_dup.py: per-library x per-tag-gene duplicate counts under the one-end
# rule, read 1 only for paired-end. See the docstring of 07_gene_dup.py. Counts only.
# Usage: sbatch --job-name=<n> --array=<tasks> --output=<exec>/logs/%A_%a_<n>_<label>.out \
#          07_gene_dup.sbatch <libraries.tsv> <regions_dir> <restricted_out_dir>
set -euo pipefail
libs=$1; regions=$2; out=$3
here=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/dup
export PYTHONNOUSERSITE=1

# columns of libraries.tsv: task cohort run individual_id bam star_log (individual_id not used)
IFS=$'\t' read -r task cohort run bam < <(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == i {print $1 "\t" $2 "\t" $3 "\t" $5}' "$libs")
[ "$task" = "$SLURM_ARRAY_TASK_ID" ] || { echo "task $SLURM_ARRAY_TASK_ID not in $libs" >&2; exit 1; }
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python "$here/07_gene_dup.py" \
  --bam "$bam" --run "$run" --cohort "$cohort" --regions "$regions" --out "$out"
