#!/bin/bash
# SS S1 mate-order probe (PRESPEC_SS amendment A1). Aligns the first 1,000,000 read pairs of each run of one library
# with the release STAR index and flags (one @RG per run, ENA mate labels as given, unsorted BAM), then counts per
# run the uniquely mapped (MAPQ 255) primary read-1 alignments inside 10 liver genes (s1_probe_genes.tsv, gene
# bodies) on the gene strand (sense) and the opposite strand. Call per run: read-1 sense fraction > 0.9 -> swap
# (ENA _2 becomes read 1), < 0.1 -> keep, else ambiguous; < 500 counted reads -> ambiguous (amendment A1.1; was 2,000).
# Usage: s1_probe.sh <S> <scratch_dir> <run>...   Writes the table (run, sense, antisense, sense_frac, call) to stdout.
set -euo pipefail
S=$1; T=$2; shift 2
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
FQ=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/data/external/SteatoSite/fastq
BIN=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
IDX=/gpfs/commons/home/jameslee/reference_genome/gencode_v49/star_index
r1=$(for r in "$@"; do printf "%s," "$FQ/${r}_1.fastq.gz"; done | sed 's/,$//')
r2=$(for r in "$@"; do printf "%s," "$FQ/${r}_2.fastq.gz"; done | sed 's/,$//')
rg=$(for r in "$@"; do printf "ID:%s SM:%s , " "$r" "$S"; done | sed 's/ , $//')
$BIN/STAR --runThreadN 8 --genomeDir $IDX --readFilesIn "$r1" "$r2" --readFilesCommand bash $SS/scripts/s1_head4m.sh \
  --outSAMattrRGline $rg --outFileNamePrefix "$T/probe." --outSAMtype BAM Unsorted > "$T/probe.star.log" 2>&1
$BIN/samtools view -q 255 -f 64 -F 2308 "$T/probe.Aligned.out.bam" \
  | awk -F'\t' -v G=$SS/scripts/s1_probe_genes.tsv -v runs="$*" '
    BEGIN { while ((getline l < G) > 0) { split(l, a, "\t"); n++; c[n]=a[1]; lo[n]=a[2]; hi[n]=a[3]; st[n]=a[4] }
            m = split(runs, R, " "); for (i = 1; i <= m; i++) { s[R[i]] = 0; x[R[i]] = 0 } }
    { rg = ""; for (i = 12; i <= NF; i++) if (substr($i, 1, 5) == "RG:Z:") rg = substr($i, 6)
      for (k = 1; k <= n; k++) if ($3 == c[k] && $4 >= lo[k] && $4 <= hi[k]) {
        rev = int($2 / 16) % 2; if ((rev == 0 && st[k] == "+") || (rev == 1 && st[k] == "-")) s[rg]++; else x[rg]++; break } }
    END { for (i = 1; i <= m; i++) { r = R[i]; t = s[r] + x[r]; f = t ? s[r] / t : -1
            call = (t < 500) ? "ambiguous" : (f > 0.9 ? "swap" : (f < 0.1 ? "keep" : "ambiguous"))
            printf "%s\t%d\t%d\t%.4f\t%s\n", r, s[r], x[r], f, call } }'
