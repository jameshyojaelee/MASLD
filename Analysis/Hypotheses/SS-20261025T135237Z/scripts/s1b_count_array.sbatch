#!/bin/bash
#SBATCH --job-name=SteatoSITE
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=io,cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-657%100
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/count/%a_%A.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z/logs/count/%a_%A.err
# SS job S1b: r8 verified fragment counts per SteatoSITE library, using the unchanged release counter
# (scripts/analysis/alt_contig_recount/run_verified_count.sbatch: -s 2 -M -R BAM -p --countReadPairs -B, r8 GTF).
# Runs only if S0a found read 1 antisense in both smoke samples. Node-local scratch is removed on exit.
set -euo pipefail
SS=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/SS-20261025T135237Z
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
grep -q "read1 antisense in both (>0.8): yes" $SS/out_s0a_smoke/strand_gate.txt || { echo "strand gate not passed"; exit 1; }
S=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $SS/scripts/s1_samples.txt)
export SAMPLE_ID=$S LIBRARY_LAYOUT=PAIRED
export SOURCE_BAM=$ROOT/RNA-seq/Human/Patient_Cohorts/results/SteatoSITE/alignments/star/$S/$S.Aligned.sortedByCoord.out.bam
export ANNOTATION_ROOT=$ROOT/figures/candidates/resource-corrections-20261005-v1/annotation_resolution_r8
export COUNT_OUTPUT_ROOT=$SS/out_s1b_counts/samples/$S
mkdir -p $SS/out_s1b_counts/samples
if [[ -s "$COUNT_OUTPUT_ROOT/count_completion.json" ]]; then echo "already counted: $S"; exit 0; fi
[ -e "$COUNT_OUTPUT_ROOT" ] && mv "$COUNT_OUTPUT_ROOT" "${COUNT_OUTPUT_ROOT}_failed_${SLURM_JOB_ID}"
unset SLURM_TMPDIR
export COUNT_SCRATCH_ROOT=$(mktemp -d "/scratch/ss-count-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}-XXXXXX")
trap 'rm -rf -- "$COUNT_SCRATCH_ROOT"' EXIT
bash $ROOT/scripts/analysis/alt_contig_recount/run_verified_count.sbatch
echo "counted $S"
