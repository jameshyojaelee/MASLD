#!/bin/bash
#SBATCH --job-name=liver_xspecies_p2_5
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/concordance_p2_5_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/logs/concordance_p2_5_%j.err
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00

set -euo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

SCRIPT_DIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Cross_Species_Concordance/scripts"

echo "=== Resuming Pipeline from Phase 2 ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo ""

# Phase 2: Pathway-Level Concordance
echo "================================================================"
echo "PHASE 2: Pathway-Level Concordance"
echo "================================================================"
echo "[Phase 2a] FGSEA concordance..."
Rscript "$SCRIPT_DIR/02a_fgsea_concordance.R"
echo ""
echo "[Phase 2b] ssGSEA concordance..."
Rscript "$SCRIPT_DIR/02b_ssgsea_concordance.R"
echo ""
echo "[Phase 2c] ORA concordance..."
Rscript "$SCRIPT_DIR/02c_ora_concordance.R"
echo ""

# Phase 3: Network-Level Concordance
echo "================================================================"
echo "PHASE 3: Network-Level Concordance"
echo "================================================================"
echo "[Phase 3a] WGCNA module preservation..."
Rscript "$SCRIPT_DIR/03a_wgcna_module_preservation.R"
echo ""
echo "[Phase 3b] TF + PROGENy activity concordance..."
Rscript "$SCRIPT_DIR/03b_tf_activity_concordance.R"
echo ""

# Phase 4: Visualizations
echo "================================================================"
echo "PHASE 4: Visualizations"
echo "================================================================"
Rscript "$SCRIPT_DIR/04_concordance_visualizations.R"
echo ""

# Phase 5: Unified Atlas
echo "================================================================"
echo "PHASE 5: Unified Concordance Atlas"
echo "================================================================"
Rscript "$SCRIPT_DIR/05_unified_concordance_atlas.R"
echo ""

echo "=== Pipeline Complete ==="
echo "Date: $(date)"
echo ""
echo "Results:"
ls -la "$SCRIPT_DIR/../results/"
echo ""
echo "Plots:"
ls -la "$SCRIPT_DIR/../plots/"
