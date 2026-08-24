#!/usr/bin/env bash
#SBATCH --job-name=masld-dnabert2-convert
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
WRAPPER="${ROOT}/slurm/dna_lm_convert_dnabert2_weights_only.sbatch"
exec bash "${ROOT}/scripts/dna_lm_convert_weights_only_driver.sh" \
  dnabert2 \
  dnabert2.pytorch_model.bin \
  7ff39ec77a484dd01070a41bfd6e95cdd7247bec80fe357ab43a4be33687aeba \
  468354983 \
  "${WRAPPER}"
