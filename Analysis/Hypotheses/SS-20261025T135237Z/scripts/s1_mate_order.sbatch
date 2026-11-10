#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/%x_%j.err
# SS S1 mate order (PRESPEC_SS amendment A1, revised). ENA read 1/read 2 labels are swapped in some runs (S1169: lanes
# of one library differ). Per run, the first 500,000 read-1 sequences are matched to exact 25-mers of 10 liver
# secretory genes (canonical exons, out_s1_mate_order_kmers/) in sense and antisense; read1_sense_frac = sense hits /
# all hits. Call: > 0.6 read 1 sense -> mates swapped for STAR (read 1 antisense, kit orientation, r8 -s 2);
# < 0.4 kept; otherwise ambiguous (library stops). Strand only; no expression value is used.
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
FQ=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq
K=$SS/out_s1_mate_order_kmers
one() { r=$1; t=$(mktemp -p "$K" .r1_XXXXXX); zcat "$FQ/${r}_1.fastq.gz" | head -2000000 | awk 'NR%4==2' > "$t" || true
  s=$(grep -c -F -f "$K/kmers_sense.txt" "$t" || true); a=$(grep -c -F -f "$K/kmers_antisense.txt" "$t" || true); rm -f "$t"
  awk -v r=$r -v s=$s -v a=$a 'BEGIN{f=(s+a)?s/(s+a):-1; c=(s+a<1000)?"ambiguous":(f>0.6?"swap":(f<0.4?"keep":"ambiguous")); printf "%s\t%d\t%d\t%.4f\t%s\n", r, s, a, f, c}'; }
export -f one; export FQ K
awk -F'\t' 'NR>1{print $1}' "$FQ/PRJEB58625_run_metadata.tsv" | xargs -P 16 -I{} bash -c 'one {}' > "$K/run_mate_order.tmp"
{ printf 'run_accession\tread1_sense_hits\tread1_antisense_hits\tread1_sense_frac\tcall\n'; sort "$K/run_mate_order.tmp"; } > "$K/run_mate_order.tsv"; rm "$K/run_mate_order.tmp"
cut -f5 "$K/run_mate_order.tsv" | tail -n +2 | sort | uniq -c
