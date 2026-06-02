#!/bin/bash
#SBATCH --job-name=masld_ds_dl_array
#SBATCH --partition=cpu
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-6
#SBATCH --output=logs/download_array_%A_%a.out
#SBATCH --error=logs/download_array_%A_%a.err

date

BASE_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/References"
EQTL_DIR="${BASE_DIR}/eQTL/MASLD_sc_eQTL_2025"
PQTL_DIR="${BASE_DIR}/pQTL/UKB_PPP_2023"

mkdir -p ${EQTL_DIR}
mkdir -p ${PQTL_DIR}
mkdir -p logs

echo "Executing SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

case ${SLURM_ARRAY_TASK_ID} in
  1)
    echo "Task 1: Downloading MASLD sc-eQTL Meta_SeuratList.rds (Zenodo)"
    cd ${EQTL_DIR}
    wget -q --show-progress "https://zenodo.org/records/14867962/files/1.Meta_SeuratList.rds?download=1" -O 1.Meta_SeuratList.rds
    ;;
  2)
    echo "Task 2: Downloading MASLD sc-eQTL Expression_Covar.zip (Zenodo)"
    cd ${EQTL_DIR}
    wget -q --show-progress "https://zenodo.org/records/14867962/files/2.Expression_Covar.zip?download=1" -O 2.Expression_Covar.zip
    unzip -o 2.Expression_Covar.zip
    rm 2.Expression_Covar.zip
    ;;
  3)
    echo "Task 3: Downloading MASLD sc-eQTL sc_eQTL_summary.zip (Zenodo)"
    cd ${EQTL_DIR}
    wget -q --show-progress "https://zenodo.org/records/14867962/files/3.sc_eQTL_summary.zip?download=1" -O 3.sc_eQTL_summary.zip
    unzip -o 3.sc_eQTL_summary.zip
    rm 3.sc_eQTL_summary.zip
    ;;
  4)
    echo "Task 4: Downloading MASLD sc-eQTL Bulk_RNA_matrix.txt (Zenodo)"
    cd ${EQTL_DIR}
    wget -q --show-progress "https://zenodo.org/records/14867962/files/4.Bulk_RNA_matrix.txt?download=1" -O 4.Bulk_RNA_matrix.txt
    ;;
  5)
    echo "Task 5: Downloading MASLD sc-eQTL hdWGCNA_module.txt (Zenodo)"
    cd ${EQTL_DIR}
    wget -q --show-progress "https://zenodo.org/records/14867962/files/5.hdWGCNA_module.txt?download=1" -O 5.hdWGCNA_module.txt
    ;;
  6)
    echo "Task 6: Downloading UKB-PPP pQTL summary statistics (AWS S3 sync)"
    cd ${PQTL_DIR}
    module load awscli || echo "awscli module not found. Assuming aws command is in PATH."
    aws s3 sync s3://ukbiobank.opendata.sagebase.org/ . --no-sign-request --exclude "*" --include "*summary_stats*" || {
        echo "WARNING: AWS sync failed. The bucket might require Synapse authentication/signed requests."
        echo "Manual download may be required via Synapse (syn51364943)."
    }
    ;;
  *)
    echo "Unknown task ID: ${SLURM_ARRAY_TASK_ID}"
    exit 1
    ;;
esac

echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID} finished."
date
