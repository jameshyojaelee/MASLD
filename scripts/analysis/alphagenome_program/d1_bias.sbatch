#!/usr/bin/env bash
#SBATCH --job-name=agp_d1_bias
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_d1_bias_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_d1_bias_%j.err
set -euo pipefail
PROJECT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=${OUT:?set OUT}
SCRIPTS=$PROJECT/scripts/analysis/alphagenome_program
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
module load SAMtools/1.23.1
module load BWA/0.7.18-GCCcore-12.3.0
module load Bowtie2/2.5.4-linux-x86_64   # the exact version named in the GSE244832 BAM @PG line
echo "host=$(hostname) start=$(date -u +%FT%TZ)"
bwa 2>&1 | grep -m1 Version || true
bowtie2 --version | head -1

"$PY" "$SCRIPTS/d1_04_mapping_bias.py" "$OUT" 16
"$PY" "$SCRIPTS/d1_05_manifest.py" "$OUT"
echo "end=$(date -u +%FT%TZ)"
