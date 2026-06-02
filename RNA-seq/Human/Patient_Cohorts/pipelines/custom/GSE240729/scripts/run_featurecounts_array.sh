#!/bin/bash
#SBATCH --job-name=fc_GSE240729
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --array=1-66
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/logs/fc_array_%A_%a.log

set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

SCRIPTDIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE240729/scripts"
OUTDIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE240729/counts/featurecounts"
GTF="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
mkdir -p "${OUTDIR}/per_sample"

BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SCRIPTDIR}/bam_list.txt")
SAMPLE=$(basename "$(dirname "$BAM")")

echo "Task ${SLURM_ARRAY_TASK_ID}: ${SAMPLE}"

featureCounts \
    -T 4 \
    -p --countReadPairs \
    -s 0 \
    -a "${GTF}" \
    -o "${OUTDIR}/per_sample/${SAMPLE}.counts.txt" \
    "${BAM}"

echo "Done: $(date)"
