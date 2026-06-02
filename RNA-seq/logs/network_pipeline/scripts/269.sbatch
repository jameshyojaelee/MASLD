#!/bin/bash
#SBATCH --job-name=net_269_export
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00


eval "$(micromamba shell hook -s bash)"
micromamba activate spatial
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
echo "269 269_export_portal_json.py | Start: $(date)"
python 269_export_portal_json.py
echo "269 269_export_portal_json.py | Exit: $? | End: $(date)"
