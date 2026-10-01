#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=90:00:00
# Job name, array range and log path are set on the sbatch command line (numbers only).
#
# B-Dup step 3, one library per array task (task = the "task" column of libraries.tsv).
# Duplicate-read fraction at tag-gene exons. Reads no allele and writes counts only.
#   1. Primary, uniquely mapped reads (STAR MAPQ 255; not unmapped, secondary or
#      supplementary) overlapping the merged tag-gene bodies -> temporary BAM. Gene
#      bodies, not exons, so both mates of a spliced pair are kept together.
#   2. samtools collate | fixmate -m | sort | markdup -c (earlier flags cleared).
#      A mate lying outside the gene bodies is absent. fixmate (samtools 1.22.1) keeps
#      0x1 on such a read but gives it no MC tag and sets RNEXT to '*'; markdup then
#      treats it as a single read. It is counted in exon_mate_absent_reads as a paired
#      read without an MC tag (checked on SRR9036307, job 21974415: 16,698 such reads
#      in the subset = markdup SINGLE; 0 reads carry 0x8).
#      Duplicate rules: complete pairs use the pair rule (both mate ends must match).
#      Single reads, which here means single-end reads and mate-absent paired reads, use
#      the one-end rule (unclipped 5' position and strand), and markdup compares them
#      with the ends of complete pairs as well as with other single reads. In paired-end
#      libraries the ~1% of exon reads that are mate-absent are therefore marked more
#      often than pairs (R-DATA review: SRR21623070 81.4% vs 59.0%) and raise
#      exon_dup_fraction slightly (0.4443 vs ~0.4387 from pairs alone). exon_dup_fraction
#      is comparable only within one layout; 08_collect_v2.py adds a one-end column
#      defined the same way for both layouts and pair-only columns.
#   3. Counts of reads overlapping the merged tag-gene exons: all, duplicate-flagged,
#      paired, mate absent. Each read counted once (samtools view -M). Overlap is
#      samtools view -L overlap: the aligned span from POS to the reference end,
#      spliced-out introns included.
# Temporary BAMs sit in <restricted>/tmp/<run>/ and are removed by name at the end;
# nothing is deleted recursively.
# Inputs : libraries.tsv, tag_gene_bodies.bed, tag_gene_exons.bed (01_regions.py)
# Outputs: <restricted>/per_library/<run>.tsv (one row of counts),
#          <restricted>/markdup_stats/<run>.txt (samtools markdup -f report for the
#          gene-body subset, including ESTIMATED_LIBRARY_SIZE)
# Usage: sbatch --job-name=<n> --array=<tasks> --output=<exec>/logs/%A_%a.out 03_dup_library.sbatch <libraries.tsv> <regions_dir> <restricted_dir>
set -euo pipefail
libs=$1; regions=$2; out=$3
st=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools
bodies=$regions/tag_gene_bodies.bed
exons=$regions/tag_gene_exons.bed

# Pass only non-empty columns to read: with IFS=tab, bash collapses empty fields, and
# individual_id is empty for the 49 libraries that failed identity calling.
IFS=$'\t' read -r task cohort run bam < <(awk -F'\t' -v i="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == i {print $1 "\t" $2 "\t" $3 "\t" $5}' "$libs")
[ "$task" = "$SLURM_ARRAY_TASK_ID" ] || { echo "task $SLURM_ARRAY_TASK_ID not in $libs" >&2; exit 1; }
mkdir -p "$out/per_library" "$out/markdup_stats" "$out/tmp"
row=$out/per_library/$run.tsv
[ -s "$row" ] && { echo "$run already done"; exit 0; }
w=$out/tmp/$run; mkdir -p "$w"

# --no-PG: without it samtools adds its own @PG line (PN:samtools) to the header it prints.
# Rows written before this change (all of full-20260929T134856Z) list "samtools,STAR"; the
# samtools entry there is the viewer itself (R-DATA review; 07 rechecks with --no-PG).
pg=$("$st" view -H --no-PG "$bam" | awk -F'\t' '$1 == "@PG" {for (i = 2; i <= NF; i++) if ($i ~ /^PN:/) print substr($i, 4)}' | sort -u | paste -sd, -)

"$st" view -@ 3 -b -M -L "$bodies" -q 255 -F 0x904 -o "$w/subset.bam" "$bam"
subset_reads=$("$st" view -c "$w/subset.bam")
pre_flagged=$("$st" view -c -f 0x400 "$w/subset.bam")

"$st" collate -@ 2 -O -u "$w/subset.bam" "$w/collate" \
  | "$st" fixmate -m -u - - \
  | "$st" sort -@ 2 -m 1500M -u -T "$w/sort" - \
  | "$st" markdup -@ 2 -c -f "$out/markdup_stats/$run.txt" - "$w/marked.bam"
"$st" index "$w/marked.bam"

count() { "$st" view -c -M -L "$exons" "$@" "$w/marked.bam"; }
exon_reads=$(count)
exon_dup=$(count -f 0x400)
exon_paired=$(count -f 0x1)
exon_mate_absent=$(count -f 0x1 -e '![MC]')

printf 'cohort\trun\tpg_programs\tsubset_reads\tpreexisting_dup_flags_in_subset\texon_reads\texon_dup_reads\texon_paired_reads\texon_mate_absent_reads\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$cohort" "$run" "${pg:-none}" "$subset_reads" "$pre_flagged" "$exon_reads" "$exon_dup" "$exon_paired" "$exon_mate_absent" > "$row.part"
mv "$row.part" "$row"

# remove this task's temporaries by name (no recursive deletion)
rm -f "$w/subset.bam" "$w/marked.bam" "$w/marked.bam.bai"
rmdir "$w" 2>/dev/null || echo "note: $w not empty after cleanup: $(ls "$w")"
echo "$run: exon_reads=$exon_reads exon_dup=$exon_dup"
