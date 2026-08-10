#!/bin/bash
#SBATCH --job-name=fC_NA5120_c0
#SBATCH --partition=cpu,io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/logs/featurecounts_s2/PRJNA512027_chunk0_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/logs/featurecounts_s2/PRJNA512027_chunk0_%j.err

echo "ERROR: PRJNA512027 is retired/noncanonical; recount is prohibited." >&2
exit 64

set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

echo "=== featureCounts chunk 0: PRJNA512027 (47/185 BAMs, PE, -s 2) ==="
echo "Started: $(date)"

BAMS=$(cat /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/scripts/featurecounts_s2/PRJNA512027_chunk0_bams.txt | tr '\n' ' ')

featureCounts \
    -T 8 \
    -p -B \
    -s 2 \
    -a "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz" \
    -o "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/PRJNA512027/counts/featurecounts/gene_counts_chunk0.txt" \
    ${BAMS}

echo "Chunk 0 done: $(date)"
