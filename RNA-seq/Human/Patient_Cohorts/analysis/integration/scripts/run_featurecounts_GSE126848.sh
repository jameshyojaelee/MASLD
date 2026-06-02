#!/bin/bash
#SBATCH --job-name=liver_fc_GSE126848
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE126848/logs/featurecounts_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE126848/logs/featurecounts_%j.err
#SBATCH --partition=cpu
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --account=nslab

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts"
GTF="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
OUTDIR="${BASE}/results/GSE126848/counts/featurecounts"
BAMS=$(find ${BASE}/results/GSE126848/alignments/star -name "*.Aligned.sortedByCoord.out.bam" | sort | tr '\n' ' ')

mkdir -p "${OUTDIR}"

echo "Running featureCounts on 57 GSE126848 BAMs..."
featureCounts \
  -T 8 \
  -s 0 \
  -a "${GTF}" \
  -o "${OUTDIR}/gene_counts.txt" \
  ${BAMS}

echo "featureCounts complete for GSE126848"
