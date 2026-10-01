#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=90:00:00
# Job name and log path are set on the sbatch command line (numbers only).
#
# B-Dup step 2: do the STAR BAMs already carry duplicate flags (0x400)?
# For the first <n_per_cohort> libraries of each cohort in libraries.tsv (sorted by run),
# runs samtools flagstat on the whole BAM and lists the @PG program names in the header.
# Reads no allele. Outputs counts only.
# Inputs : libraries.tsv from 01_regions.py
# Outputs: <restricted>/flag_check_per_library.tsv (cohort, run, pg_programs, total,
#          primary, secondary, supplementary, duplicates, primary_duplicates,
#          paired_in_sequencing); <exec>/flag_check_per_cohort.tsv
# Usage: sbatch --job-name=<n> --output=<exec>/logs/%j.out 02_flag_check.sbatch <libraries.tsv> <n_per_cohort> <restricted_dir> <exec_dir>
set -euo pipefail
libs=$1; n=$2; restricted=$3; exec_dir=$4
st=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools

per_lib=$restricted/flag_check_per_library.tsv
[ -e "$per_lib" ] && { echo "refusing to overwrite $per_lib" >&2; exit 1; }
printf 'cohort\trun\tpg_programs\ttotal\tprimary\tsecondary\tsupplementary\tduplicates\tprimary_duplicates\tpaired_in_sequencing\n' > "$per_lib.part"

# columns of libraries.tsv: task cohort run individual_id bam star_log
awk -F'\t' -v n="$n" 'NR > 1 && ++seen[$2] <= n {print $2 "\t" $3 "\t" $5}' "$libs" |
while IFS=$'\t' read -r cohort run bam; do
  # --no-PG: otherwise samtools adds its own @PG line; the 2026-09-27 run (job 21941401) did
  # not use it, so its "samtools" entry is the viewer itself, not a program run on the BAM.
  pg=$("$st" view -H --no-PG "$bam" | awk -F'\t' '$1 == "@PG" {for (i = 2; i <= NF; i++) if ($i ~ /^PN:/) print substr($i, 4)}' | sort -u | paste -sd, -)
  "$st" flagstat -@ 8 -O tsv "$bam" > "$restricted/flag_check.$run.tmp"
  # flagstat tsv rows: QC-passed, QC-failed, label; report passed + failed
  get() { awk -F'\t' -v k="$1" '$3 == k {print $1 + $2}' "$restricted/flag_check.$run.tmp"; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$cohort" "$run" "${pg:-none}" \
    "$(get 'total (QC-passed reads + QC-failed reads)')" "$(get primary)" "$(get secondary)" \
    "$(get supplementary)" "$(get duplicates)" "$(get 'primary duplicates')" "$(get 'paired in sequencing')" >> "$per_lib.part"
  rm -f "$restricted/flag_check.$run.tmp"
  echo "$cohort $run done"
done
mv "$per_lib.part" "$per_lib"

awk -F'\t' 'NR > 1 {
    n[$1]++; dup[$1] += $8; tot[$1] += $4; paired[$1] += ($10 > 0)
    if ($3 ~ /[Mm]ark[Dd]up|MarkDuplicates|[Dd]edup|UMI/) pgdup[$1]++
    progs[$1] = $3
  }
  END {
    print "cohort\tn_libraries_checked\ttotal_records\tduplicate_flagged_records\tn_paired_libraries\tn_headers_naming_a_dup_tool\tpg_programs_example"
    for (c in n) print c "\t" n[c] "\t" tot[c] "\t" dup[c] "\t" paired[c] "\t" pgdup[c] + 0 "\t" progs[c]
  }' "$per_lib" | (read -r h; echo "$h"; sort) > "$exec_dir/flag_check_per_cohort.tsv"
"$st" --version | head -2 > "$exec_dir/flag_check_samtools_version.txt"
cat "$exec_dir/flag_check_per_cohort.tsv"
