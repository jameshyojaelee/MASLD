#!/bin/bash
#SBATCH --job-name=gencodeid
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/noncoding_resource/slurm_logs/gencodeid_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/noncoding_resource/slurm_logs/gencodeid_%j.err

set -euo pipefail

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CANDIDATE_ROOT="${PROJECT_ROOT}/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11"
GTF="${PROJECT_ROOT}/RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/inputs/BG001-DECISION/external_reference/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
OUTPUT="${CANDIDATE_ROOT}/annotation/GENCODE_v49/run_${SLURM_JOB_ID}"

python3 "${PROJECT_ROOT}/scripts/manuscript/noncoding_resource/build_gene_identity.py" \
    --gtf "${GTF}" \
    --output-dir "${OUTPUT}"
