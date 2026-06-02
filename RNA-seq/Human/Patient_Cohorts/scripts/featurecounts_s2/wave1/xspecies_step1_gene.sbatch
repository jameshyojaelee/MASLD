#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=fgsea
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/xspecies_step1_gene_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/wave1/xspecies_step1_gene_%j.err
# Wave1 [E] step1: 01b (nafl-specific sig from -s2 merged_counts_raw) -> 01 (4x5 gene concordance,
# reads -s2 dream_results.csv + mouse per-diet DE). Serial: 01 needs 01b output.
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
S=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "=== xspecies 01b generate nafl-specific | $(date) ==="
Rscript "$S/01b_generate_nafl_specific.R"
if [ $? -ne 0 ]; then echo "FAILED: 01b_generate_nafl_specific"; exit 1; fi
echo "=== xspecies 01 corrected gene concordance | $(date) ==="
Rscript "$S/01_corrected_gene_concordance.R"
if [ $? -ne 0 ]; then echo "FAILED: 01_corrected_gene_concordance"; exit 1; fi
echo "=== done $(date) ==="
