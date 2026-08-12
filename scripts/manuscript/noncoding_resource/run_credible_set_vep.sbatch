#!/bin/bash
#SBATCH --job-name=ncvep
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncvep_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/ncvep_%j.err

set -euo pipefail
umask 027

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
MEMBERS="${PROJECT_ROOT}/GWAS/finemapping/results/candidates/noncoding-dna-precoloc-2026-08-12/run_19780781/credible_set_members_annotated.tsv"
OUTPUT="${NONCODING_VEP_OUTPUT:-${PROJECT_ROOT}/GWAS/finemapping/results/candidates/noncoding-dna-vep-2026-08-12/run_${SLURM_JOB_ID}}"

test ! -e "${OUTPUT}"

python3 scripts/manuscript/noncoding_resource/annotate_credible_set_vep.py \
  --members "${MEMBERS}" \
  --output-dir "${OUTPUT}" \
  --batch-size 200
