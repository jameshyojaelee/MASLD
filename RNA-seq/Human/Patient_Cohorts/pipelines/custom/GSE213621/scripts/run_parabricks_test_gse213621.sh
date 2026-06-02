#!/bin/bash
#SBATCH --job-name=liver_parabricks
#SBATCH --partition=gpu
#SBATCH --gres=gpu:L40S:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE213621/slurm_logs/pb_test_GSE213621_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE213621/slurm_logs/pb_test_GSE213621_%j.err

# Load Apptainer for the Parabricks wrapper
module load apptainer/1.2.3

# Ensure pbrun wrapper is in PATH
export PATH="/gpfs/commons/home/jameslee/.local/bin:$PATH"

# Paths
REFERENCE="/gpfs/commons/home/jameslee/reference_genome/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
GTF="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
STAR_INDEX="/gpfs/commons/home/jameslee/reference_genome/gencode_v49/star_index"

OUTDIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/results/GSE213621/parabricks_run"

# Sample 1
SAMPLE1="SRR21622825"
FQ1_1="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE213621/fastq/SRR21622825_1.fastq.gz"
FQ1_2="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE213621/fastq/SRR21622825_2.fastq.gz"

# Sample 2
SAMPLE2="SRR21622826"
FQ2_1="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE213621/fastq/SRR21622826_1.fastq.gz"
FQ2_2="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/data/raw/GSE213621/fastq/SRR21622826_2.fastq.gz"

echo "Creating output directories..."
mkdir -p ${OUTDIR}/${SAMPLE1}
mkdir -p ${OUTDIR}/${SAMPLE2}

echo "Running Parabricks rna_fq2bam for $SAMPLE1..."
pbrun rna_fq2bam \
  --ref ${REFERENCE} \
  --genome-lib-dir ${STAR_INDEX} \
  --in-fq ${FQ1_1} ${FQ1_2} \
  --read-files-command zcat \
  --num-threads 8 \
  --memory-limit 80 \
  --output-dir ${OUTDIR}/${SAMPLE1} \
  --tmp-dir ${OUTDIR}/${SAMPLE1} \
  --out-bam ${OUTDIR}/${SAMPLE1}/${SAMPLE1}.Aligned.sortedByCoord.out.bam \
  --read-group-sm ${SAMPLE1}

echo "Running Parabricks rna_fq2bam for $SAMPLE2..."
pbrun rna_fq2bam \
  --ref ${REFERENCE} \
  --genome-lib-dir ${STAR_INDEX} \
  --in-fq ${FQ2_1} ${FQ2_2} \
  --read-files-command zcat \
  --num-threads 8 \
  --memory-limit 80 \
  --output-dir ${OUTDIR}/${SAMPLE2} \
  --tmp-dir ${OUTDIR}/${SAMPLE2} \
  --out-bam ${OUTDIR}/${SAMPLE2}/${SAMPLE2}.Aligned.sortedByCoord.out.bam \
  --read-group-sm ${SAMPLE2}

echo "Parabricks testing complete."
