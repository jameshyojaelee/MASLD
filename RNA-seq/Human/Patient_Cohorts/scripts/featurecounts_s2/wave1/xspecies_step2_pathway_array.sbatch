#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-2
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/xspecies_step2_pathway_%a_%A.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/xspecies_step2_pathway_%a_%A.err
# Wave1 [E] step2: 02a fgsea / 02b ssgsea / 02c ora concordance (independent -> array 0-2).
# All read -s2 human DE + mouse DE (02b reads mouse corrected_logcpm).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
SCRIPTS=(02a_fgsea_concordance.R 02b_ssgsea_concordance.R 02c_ora_concordance.R)
THIS=${SCRIPTS[$SLURM_ARRAY_TASK_ID]}
S=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== xspecies step2 task ${SLURM_ARRAY_TASK_ID}: ${THIS} | $(date) ==="
Rscript "$S/${THIS}"
if [ $? -ne 0 ]; then echo "FAILED: ${THIS}"; exit 1; fi
echo "=== done ${THIS} $(date) ==="
