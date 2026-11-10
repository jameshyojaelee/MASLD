#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=3-00:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s1d/test_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/s1d/test_%j.err
# SS S1d test: s1d_ig_one.sh on release BAM GSE130970 SRR9036307; compare ig_fragments.tsv with the stored full-pass
# output of extract_reads_one_library.py (READS-20261009T141017Z); time and peak memory from /usr/bin/time -v.
set -euo pipefail
V=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/verify/s1d-20261025T173807Z
/usr/bin/time -v /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s1d_ig_one.sh SRR9036307 /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results/GSE130970/alignments/star/SRR9036307/SRR9036307.Aligned.sortedByCoord.out.bam $V/ig 2> $V.time || { cat $V.time; exit 1; }
cat $V.time | grep -E "Elapsed|Maximum resident"
REF=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/READS-20261009T141017Z/libraries/GSE130970/SRR9036307/ig_fragments.tsv
cmp $V/ig/SRR9036307/ig_fragments.tsv $REF && echo "ig_fragments.tsv IDENTICAL to stored full pass"
cat $V/ig/SRR9036307/library_qc.tsv
/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s1d_ig_one.sh SRR9036307 /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/RNA-seq/Human/Patient_Cohorts/results/GSE130970/alignments/star/SRR9036307/SRR9036307.Aligned.sortedByCoord.out.bam $V/ig && echo "OVERWRITE NOT REFUSED" || echo "rerun refused as intended"
