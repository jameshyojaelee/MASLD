#!/bin/bash
#SBATCH --job-name=s3d_loco_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3d_loco_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/SS-20261025T135237Z/logs/s3d_loco_20261118_%j.err

set -euo pipefail
PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
LOCO_DIR="${PROJECT_ROOT}/Analysis/Deconvolution/candidates/loco_20261118"
OUT_S3D="${PROJECT_ROOT}/Analysis/Hypotheses/SS-20261025T135237Z/out_s3d_loco_20261118"
rm -rf "${OUT_S3D}"

/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla \
  "${PROJECT_ROOT}/Analysis/Hypotheses/SS-20261025T135237Z/scripts/s3d_programs_20261118.R" \
  --input "${PROJECT_ROOT}/Analysis/Hypotheses/SS-20261025T135237Z/out_s2_qc" \
  --stage "${PROJECT_ROOT}/data/external/SteatoSite/fastq/PRJEB58625_fibrosis_stage_20261025T133450Z.tsv" \
  --out "${OUT_S3D}" \
  --composition "${LOCO_DIR}/standard_steatosite_ilr8.tsv" \
  --loco_dir "${LOCO_DIR}/chr" \
  --h54 skip
