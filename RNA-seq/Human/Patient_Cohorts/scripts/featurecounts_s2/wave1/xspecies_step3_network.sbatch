#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/xspecies_step3_network_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/xspecies_step3_network_%j.err
# Wave1 [E] step3: 03a WGCNA module preservation + 03b TF/PROGENy concordance (after step1).
# Read -s2 merged_dge.rds + dream_results.csv + mouse corrected_logcpm. Serial (both write to RES).
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
S=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== xspecies 03a WGCNA module preservation | $(date) ==="
Rscript "$S/03a_wgcna_module_preservation.R"
if [ $? -ne 0 ]; then echo "FAILED: 03a_wgcna_module_preservation"; exit 1; fi
echo "=== xspecies 03b TF activity concordance | $(date) ==="
Rscript "$S/03b_tf_activity_concordance.R"
if [ $? -ne 0 ]; then echo "FAILED: 03b_tf_activity_concordance"; exit 1; fi
echo "=== done $(date) ==="
