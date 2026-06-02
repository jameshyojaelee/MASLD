#!/bin/bash
#SBATCH --job-name=liver_launcher_GSE167523
#SBATCH --output=launcher_%j.out
#SBATCH --error=launcher_%j.err
#SBATCH --partition=cpu
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

echo "Pipeline launcher started."
# Ensure MICROMAMBA is available
MICROMAMBA="/gpfs/commons/home/jameslee/.local/bin/micromamba"
if [ -x "$MICROMAMBA" ]; then
  eval "$($MICROMAMBA shell hook --shell bash)"
  micromamba activate rnaseq
else
  echo "Error: Micromamba not found at $MICROMAMBA"
  exit 1
fi

# Run the submit script
echo "Running submit_pipeline.sh..."
bash submit_pipeline.sh
