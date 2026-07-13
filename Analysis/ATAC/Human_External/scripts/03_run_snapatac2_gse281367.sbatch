#!/bin/bash
# 03_run_snapatac2_gse281367.sbatch  --  Layer A A3 (GSE281367 SnapATAC2 processing)
# Reuse the tested GSE244832 processing script (02_snapatac2_processing.py) on the
# new cohort: import fragments -> QC (min_tsse=4.0, min_frags=1000) -> tile matrix ->
# spectral -> harmony(donor) -> leiden -> marker annotation. Skip per-cell DA.
#
#SBATCH --job-name=snapatac2
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=Analysis/ATAC/Human_External/logs/snap_%j.out
#SBATCH --error=Analysis/ATAC/Human_External/logs/snap_%j.err
set -eo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate snapatac2

echo "[$(date)] SnapATAC2 processing GSE281367 on $(hostname)"
python Analysis/ATAC/Human_Multiome/scripts/02_snapatac2_processing.py \
  --input_dir  "$PROJ/Analysis/ATAC/Human_External/frag_links" \
  --output_dir "$PROJ/Analysis/ATAC/Human_External/snapatac2" \
  --sample_csv "$PROJ/data/GSE281367/metadata/donor_pairing.csv" \
  --n_workers 16 --min_tsse 4.0 --min_frags 1000 --skip_da

echo "[$(date)] Done. Outputs:"
ls -la "$PROJ/Analysis/ATAC/Human_External/snapatac2/" 2>/dev/null
