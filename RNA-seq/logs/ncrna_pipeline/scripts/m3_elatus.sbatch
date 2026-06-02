#!/bin/bash
#SBATCH --job-name=ncRNA_M3_elatus
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00

eval "$(micromamba shell hook -s bash)"
micromamba activate rapids_singlecell
export LD_LIBRARY_PATH=$(python -c 'import nvidia.cusparselt.lib; print(nvidia.cusparselt.lib.__path__[0])' 2>/dev/null):${LD_LIBRARY_PATH:-}
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
echo "Module 3: ELATUS/scRNA | Start: $(date)"
python Analysis/SingleCell/scripts/elatus_lncrna_profiling.py
echo "Module 3: Exit code: $? | End: $(date)"
