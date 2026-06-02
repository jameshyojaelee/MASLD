#!/bin/bash
#SBATCH --job-name=A4_oliva_sbeqtl
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/logs/A4_oliva_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/logs/A4_oliva_%j.err

set -eo pipefail

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
cd "$PROJ"

# Activate rnaseq env (binary lives in ~/.local/bin/micromamba; profile.d dir only has mamba.sh).
# rnaseq env's binutils activate hook references unbound ADDR2LINE; relax `-u` around activate.
set +u
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

mkdir -p "$PROJ/RNA-seq/scripts/logs"

echo "===== A4: Oliva 2020 sb-eQTL integration ====="
echo "Host : $(hostname)"
echo "Date : $(date)"
echo "Job  : ${SLURM_JOB_ID:-local}"
echo "CPUs : ${SLURM_CPUS_PER_TASK:-8}"
echo

# Re-fetch source files if missing (idempotent)
SB_DIR="$PROJ/data/gtex_sb_eqtl"
mkdir -p "$SB_DIR"

if [[ ! -f "$SB_DIR/GTEx_Analysis_v8_sbeQTLs/GTEx_Analysis_v8_sbeQTLs.txt" ]]; then
  echo "[fetch] sb-eQTL tarball missing, downloading ..."
  curl -sL --fail --max-time 600 \
    -o "$SB_DIR/GTEx_Analysis_v8_sbeQTLs.tar.gz" \
    "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/sb-eqtl/GTEx_Analysis_v8_sbeQTLs.tar.gz"
  tar -xzf "$SB_DIR/GTEx_Analysis_v8_sbeQTLs.tar.gz" -C "$SB_DIR"
fi

if [[ ! -f "$SB_DIR/GTEx_Analysis_v8_sbgenes/signif.sbgenes.txt" ]]; then
  echo "[fetch] sb-genes tarball missing, downloading ..."
  curl -sL --fail --max-time 600 \
    -o "$SB_DIR/GTEx_Analysis_v8_sbgenes.tar.gz" \
    "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/sb-eqtl/GTEx_Analysis_v8_sbgenes.tar.gz"
  tar -xzf "$SB_DIR/GTEx_Analysis_v8_sbgenes.tar.gz" -C "$SB_DIR"
fi

echo "[run] Harmonization + integration ..."
Rscript "$PROJ/RNA-seq/scripts/oliva_sbeqtl_integration.R"

echo
echo "===== Done $(date) ====="
