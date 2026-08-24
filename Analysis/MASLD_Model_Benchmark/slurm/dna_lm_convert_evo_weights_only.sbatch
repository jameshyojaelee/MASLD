#!/usr/bin/env bash
#SBATCH --job-name=masld-evo-convert
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --qos=nslab
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark/executions/%x-%j.err

ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark
WRAPPER="${ROOT}/slurm/dna_lm_convert_evo_weights_only.sbatch"
exec bash "${ROOT}/scripts/dna_lm_convert_weights_only_driver.sh" \
  evo \
  evo.evo2_7b.pt \
  c66645929dc1b9c631f5be656da8726f38946315dc9167000a615dd626fcecf4 \
  13766621200 \
  "${WRAPPER}"
