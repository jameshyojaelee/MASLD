#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
# Job name and log path are set on the sbatch command line (numbers only).
#
# B-Dup diagnostic for one paired-end library: where are the reads whose mate is absent
# from the gene-body subset, and which flags do they carry after samtools fixmate?
# Rebuilds the 03_dup_library.sbatch subset, runs collate | fixmate -m | sort, then
# prints counts and the first alignment columns (QNAME..TLEN, fields 1-9; never SEQ)
# of mate-absent reads. Temporaries are removed by name.
# Usage: sbatch ... diag_orphan_flags.sbatch <bam> <regions_dir> <work_dir>
set -euo pipefail
bam=$1; regions=$2; w=$3
st=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/samtools
mkdir -p "$w"
"$st" view -@ 3 -b -M -L "$regions/tag_gene_bodies.bed" -q 255 -F 0x904 -o "$w/subset.bam" "$bam"
"$st" collate -@ 2 -O -u "$w/subset.bam" "$w/collate" | "$st" fixmate -m -u - - \
  | "$st" sort -@ 2 -m 1500M -T "$w/sort" -o "$w/fixed.bam" -
"$st" index "$w/fixed.bam"
ex=$regions/tag_gene_exons.bed
echo "subset reads            $("$st" view -c "$w/fixed.bam")"
echo "0x8 (mate unmapped)     $("$st" view -c -f 0x8 "$w/fixed.bam")"
echo "paired, no MC tag  ![MC]           $("$st" view -c -f 0x1 -e '![MC]' "$w/fixed.bam")"
echo "paired, no MC tag  !exists([MC])   $("$st" view -c -f 0x1 -e '!exists([MC])' "$w/fixed.bam" 2>&1 | tail -1)"
echo "paired, with MC tag [MC]           $("$st" view -c -f 0x1 -e '[MC]' "$w/fixed.bam")"
echo "no MC at exons, -M                 $("$st" view -c -M -L "$ex" -f 0x1 -e '![MC]' "$w/fixed.bam")"
echo "all at exons, -M                   $("$st" view -c -M -L "$ex" "$w/fixed.bam")"
set +o pipefail  # head closes pipes early below
echo "flags of no-MC reads:"; "$st" view -f 0x1 -e '![MC]' "$w/fixed.bam" | cut -f2 | sort | uniq -c | sort -k1,1nr | head
echo "first no-MC reads (fields 1-9):"; "$st" view -f 0x1 -e '![MC]' "$w/fixed.bam" | cut -f1-9 | head -8
echo "same names in the original BAM near those positions (fields 1-9):"
"$st" view -f 0x1 -e '![MC]' "$w/fixed.bam" | head -8 | cut -f1 > "$w/names.txt"
"$st" view -f 0x1 -e '![MC]' "$w/fixed.bam" | head -8 | awk '{print $3 ":" $4}' > "$w/regions.txt"
while read -r reg; do
  c=${reg%%:*}; p=${reg#*:}
  "$st" view "$bam" "$c:$((p > 5000 ? p - 5000 : 1))-$((p + 5000))" | grep -F -f "$w/names.txt" | cut -f1-9 || true
done < "$w/regions.txt" | sort -u
rm -f "$w/subset.bam" "$w/fixed.bam" "$w/fixed.bam.bai" "$w/names.txt" "$w/regions.txt"
rmdir "$w" 2>/dev/null || ls "$w"
