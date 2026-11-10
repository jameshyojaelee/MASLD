#!/bin/bash
#SBATCH --job-name=unique-capture
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/NUCLEUS-FIX-20261107/recount_%j.log
set -euo pipefail
if [[ "$1" == SRX* ]]; then
    python3 Analysis/Hypotheses/NUCLEUS-FIX-20261107/scripts/recount_capture.py "$1"
else
    python3 Analysis/Hypotheses/NUCLEUS-FIX-20261107/scripts/recount_gse185477.py "$1"
fi
