#!/bin/bash
#SBATCH --job-name=liver_fix_script03
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fix_script03_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/fix_script03_%j.err

# Fix for: Script 03 was assigning group_binary by position (not by sample_id match),
# corrupting 278/965 sample labels in merged_dge.rds.
# This re-runs 03 (rebuild DGE), 05 (dream), 07 (consensus DEGs), 08 (pathway analysis).

set -euo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
set +u; micromamba activate rnaseq; set -u

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts
SCRIPTS=$BASE/analysis/integration/scripts

echo "=== Script 03: Rebuild merged_dge.rds (label order fix) ==="
echo "Start: $(date)"
Rscript $SCRIPTS/03_integrate_counts.R 2>&1
echo "Script 03 complete: $(date)"

echo "=== Script 05: Dream mega-analysis (6 cohorts) ==="
Rscript $SCRIPTS/05_dream_mega_analysis.R 2>&1
echo "Script 05 complete: $(date)"

echo "=== Script 07: Consensus DEGs ==="
Rscript $SCRIPTS/07_consensus_degs.R 2>&1
echo "Script 07 complete: $(date)"

echo "=== Script 08: Pathway analysis ==="
Rscript $SCRIPTS/08_pathway_analysis.R 2>&1
echo "Script 08 complete: $(date)"

echo "=== All done: $(date) ==="
