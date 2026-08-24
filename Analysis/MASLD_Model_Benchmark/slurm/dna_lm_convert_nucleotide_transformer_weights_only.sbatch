#!/usr/bin/env bash
#SBATCH --job-name=masld-nt-convert
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
WRAPPER="${ROOT}/slurm/dna_lm_convert_nucleotide_transformer_weights_only.sbatch"
exec bash "${ROOT}/scripts/dna_lm_convert_weights_only_driver.sh" \
  nucleotide_transformer \
  nucleotide_transformer.pytorch_model.bin \
  971cb721bd8cc8134d3665b5d94b195fffb18a1480ad5925d6dcc9290bf1c17f \
  1942934337 \
  "${WRAPPER}"
