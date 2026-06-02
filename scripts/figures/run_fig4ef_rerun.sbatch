#!/bin/bash
#SBATCH --job-name=fig4ef_fix
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/fig4ef_fix_%j.log
#SBATCH --error=logs/fig4ef_fix_%j.log

set -euo pipefail
BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SCRIPTS=$BASE/scripts/figures

micromamba run -n rnaseq Rscript "$SCRIPTS/fig4e_serpine1_circuit.R"
echo "[OK] 4E"

micromamba run -n rnaseq Rscript "$SCRIPTS/fig4f_cyp3a4_zonation.R"
echo "[OK] 4F"

micromamba run -n spatial python "$SCRIPTS/fig4b_spatial_overview.py"
echo "[OK] 4B"

echo "Done. PDFs:"
ls -lh "$BASE/figures/main/fig4_validation/panels/"*.pdf 2>/dev/null
