#!/usr/bin/env bash
#SBATCH --job-name=1592
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1592_raw_layout_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/1592_raw_layout_%j.err
set -euo pipefail
umask 077
ulimit -c 0
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
module load SRA-Toolkit/3.2.1
module load cellranger-arc/2.0.2
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
exec /usr/bin/python3 "$TASK_ROOT/scripts/analysis/in_biopsy_allelic/rec_codex/codex_gse223843_read_layout_pilot.py" \
    --expected-script-sha256 53771e445f3457e134a5727ded259f83e2f44ffb987c84c5d2f3f0c7e31a100f
