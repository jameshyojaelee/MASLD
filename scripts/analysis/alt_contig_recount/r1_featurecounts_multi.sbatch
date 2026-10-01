#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --array=0-4
# R1: rerun featureCounts for one F_five cohort with exactly the flags recorded in
# its effective count file (BG001 F_five provenance) plus "-M --fraction", on the
# same BAMs in the same column order. Read-only on existing results.
# Usage: sbatch --job-name=<n> --output=<out>/logs/%A_%a.out r1_featurecounts_multi.sbatch <out_dir>
set -euo pipefail
out=$1
read -r -a cohorts <<< "${COHORTS:-GSE126848 GSE130970 GSE135251 GSE162694 GSE213621}"
c=${cohorts[$SLURM_ARRAY_TASK_ID]}
prov=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/results/remediation/bg001/bg001-fragment-v211-gencode49-20260807T195243Z/arms/F_five/provenance/effective_count_sources.tsv
fc=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/featureCounts
mkdir -p "$out/counts"

src=$(awk -F'\t' -v c="$c" '$1 == c {print $2}' "$prov")
test -s "$src"
# The first header line records the original command: "featureCounts" "arg" ... "-o" "<out>" "<bam>" ...
header=$(awk 'NR == 1' "$src")
args=$(printf '%s\n' "$header" | grep -o '"[^"]*"' | tr -d '"' | awk '
  NR == 1 {next}                                  # program name
  skip {skip = 0; next}
  $0 == "-o" || $0 == "-a" || $0 == "-T" {skip = 1; next}
  /\.bam$/ {next}
  {print}' | tr '\n' ' ')
gtf=$(printf '%s\n' "$header" | grep -o '"-a" "[^"]*"' | awk -F'"' '{print $4}')
awk 'NR == 2' "$src" | tr '\t' '\n' | awk 'NR > 6' > "$out/counts/$c.bams.txt"
echo "$c original args: $args" | tee "$out/counts/$c.args.txt"
echo "gtf: $gtf" >> "$out/counts/$c.args.txt"
test -s "$out/counts/$c.bams.txt"

# shellcheck disable=SC2046,SC2086
"$fc" -T 12 $args -M --fraction -a "$gtf" -o "$out/counts/$c.multi.txt" $(cat "$out/counts/$c.bams.txt")
"$fc" -v > "$out/counts/$c.featurecounts_version.txt" 2>&1
echo "$c done"
