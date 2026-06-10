#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/18_06i_shim_%j.out
#SBATCH --error=logs/18_06i_shim_%j.err

# ---------------------------------------------------------------------------
# Regenerate the v6 cross-pillar consensus (18) from the current pillar
# inputs, then overlay it onto the LVQW-rebuilt sex_deg_classification_v3.csv
# via the 06i shim. The Jun-8 LVQW migration re-ran modules 04/05/06 (which
# regenerated v3.csv) but did NOT re-run 18 or 06i, so v3.csv lost the
# class_v6_consensus columns the sex figure requires. This restores them.
# Honest note: the v6 statistical pillars (P2/P4/P6/P8) remain the May-29
# dream-based artifacts; under the LVQW-collapsed reality the consensus is
# expected to be degenerate (0 Strong / 0 Moderate). We run it anyway so the
# figure reflects the honest collapsed result.
# ---------------------------------------------------------------------------
set -eo pipefail

SCRIPT_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
mkdir -p "${SCRIPT_DIR}/logs"
cd "${SCRIPT_DIR}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u
micromamba activate rnaseq
set -u

echo "=== 18_consensus_v6 start: $(date) on $(hostname) ==="
Rscript 18_consensus_v6.R
echo "--- 06i v6->v3 shim ---"
Rscript 06i_v6_to_v3_shim.R
echo "=== 18 + 06i end: $(date) ==="
