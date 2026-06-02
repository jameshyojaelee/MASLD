#!/bin/bash
#SBATCH --job-name=fc_GSE174478
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/pipelines/custom/GSE174478/logs/featurecounts_%j.log

set -euo pipefail

eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq

RESULTS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE174478"
GTF="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
OUTDIR="${RESULTS}/counts/featurecounts"
mkdir -p "${OUTDIR}"

# Collect all BAMs that exist (skip SRR14551000 which failed STAR)
BAMS=()
for d in ${RESULTS}/alignments/star/*/; do
    sample=$(basename "$d")
    bam="${d}/${sample}.Aligned.sortedByCoord.out.bam"
    if [ -f "$bam" ]; then
        BAMS+=("$bam")
    fi
done

echo "Running featureCounts on ${#BAMS[@]} BAMs (paired-end)"
echo "Output: ${OUTDIR}/gene_counts.txt"

featureCounts \
    -T 8 \
    -p --countReadPairs \
    -s 0 \
    -a "${GTF}" \
    -o "${OUTDIR}/gene_counts.txt" \
    "${BAMS[@]}"

echo "Done: $(date)"
echo "Samples counted: ${#BAMS[@]}"
wc -l "${OUTDIR}/gene_counts.txt"
