#!/usr/bin/env bash
#SBATCH --job-name=agp_d1_counts
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_d1_counts_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_d1_counts_%j.err
set -euo pipefail
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${OUT:?set OUT}
SCRIPTS=$PROJECT/scripts/analysis/alphagenome_program
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
module load BCFtools/1.19-GCC-13.2.0
module load SAMtools/1.23.1
echo "host=$(hostname) start=$(date -u +%FT%TZ)"
bcftools --version | head -1; samtools --version | head -1

bash "$SCRIPTS/d1_02_extract_and_count.sh" "$OUT" 16
"$PY" "$SCRIPTS/d1_03_assemble_counts.py" "$OUT"

"$PY" -m pip freeze > "$OUT/pip_freeze_rnaseq.txt" 2>/dev/null || true
echo "end=$(date -u +%FT%TZ)"
