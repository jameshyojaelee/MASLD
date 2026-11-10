#!/bin/bash
# SSX-B intronic fragment count for one SteatoSITE library (PRESPEC_SSX_B.md): featureCounts -s 2 -p --countReadPairs -B
# on the intronic GTF (unique alignments only; default overlap rules). Usage: ssx_b_count_one.sh <sample> <out dir>
set -euo pipefail
S=$1; O=$2; mkdir -p $O
[[ -s $O/$S.intron.tsv.summary ]] && { echo "already counted: $S"; exit 0; }
BAM=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results/SteatoSITE/alignments/star/$S/$S.Aligned.sortedByCoord.out.bam
T=$(mktemp -d "/scratch/ssxb-${SLURM_JOB_ID:-x}-$S-XXXXXX"); trap 'rm -rf -- "$T"' EXIT
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/featureCounts -T ${SLURM_CPUS_PER_TASK:-4} -s 2 -p --countReadPairs -B \
  -t intron -g gene_id -a /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SSX-20261102T133212Z/out_b_gtf-20261102T133544Z/intronic_features.gtf --tmpDir $T -o $T/$S.intron.tsv $BAM
cut -f1,7 $T/$S.intron.tsv | tail -n +3 | gzip > $T/$S.intron.counts.gz
mv $T/$S.intron.counts.gz $O/ && mv $T/$S.intron.tsv.summary $O/
echo "counted $S"
