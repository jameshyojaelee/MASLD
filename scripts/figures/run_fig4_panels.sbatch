#!/bin/bash
#SBATCH --job-name=fig4_panels
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/fig4_panels_%j.log
#SBATCH --error=logs/fig4_panels_%j.log

set -euo pipefail

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=$BASE/scripts/figures
LOGDIR=$SCRIPTS/../logs
mkdir -p "$LOGDIR"

echo "=== Fig 4 panel generation ==="
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"

# ── R panels (rnaseq env) ─────────────────────────────────────────────────────
micromamba run -n rnaseq Rscript "$SCRIPTS/fig4a_proteomics_overview.R"
echo "[OK] 4A proteomics overview"

micromamba run -n rnaseq Rscript "$SCRIPTS/fig4c_nr_triad.R"
echo "[OK] 4C NR triad"

micromamba run -n rnaseq Rscript "$SCRIPTS/fig4d_hkdc1_circuit.R"
echo "[OK] 4D HKDC1 circuit"

micromamba run -n rnaseq Rscript "$SCRIPTS/fig4e_serpine1_circuit.R"
echo "[OK] 4E SERPINE1 circuit"

micromamba run -n rnaseq Rscript "$SCRIPTS/fig4f_cyp3a4_zonation.R"
echo "[OK] 4F CYP3A4 zonation"

# ── Python panel (spatial env) ────────────────────────────────────────────────
micromamba run -n spatial python "$SCRIPTS/fig4b_spatial_overview.py"
echo "[OK] 4B spatial overview"

echo ""
echo "=== All panels complete ==="
echo "Outputs: $BASE/figures/main/fig4_validation/panels/"
ls -lh "$BASE/figures/main/fig4_validation/panels/"*.pdf 2>/dev/null || echo "(no PDFs yet)"
echo "Finished: $(date)"
