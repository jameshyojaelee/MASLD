#!/bin/bash
#SBATCH --job-name=fC_213621_c2
#SBATCH --partition=cpu,io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/logs/featurecounts_s2/GSE213621_chunk2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/logs/featurecounts_s2/GSE213621_chunk2_%j.err

set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

echo "=== featureCounts chunk 2: GSE213621 (92/367 BAMs, PE, -s 2) ==="
echo "Started: $(date)"

BAMS=$(cat /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/GSE213621_chunk2_bams.txt | tr '\n' ' ')

featureCounts \
    -T 8 \
    -p -B \
    -s 2 \
    -a "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz" \
    -o "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE213621/counts/featurecounts/gene_counts_chunk2.txt" \
    ${BAMS}

echo "Chunk 2 done: $(date)"
