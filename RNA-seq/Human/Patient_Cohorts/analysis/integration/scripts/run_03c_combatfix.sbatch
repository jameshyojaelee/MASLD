#!/bin/bash
#SBATCH --job-name=combatfix
#SBATCH --partition=cpu
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/combatfix_%j.out
#SBATCH --error=logs/combatfix_%j.err
# ---------------------------------------------------------------------------
# Audit P2#20 remediation (2026-06-29): rebuild the human batch-corrected logCPM
# under the include_in_mega gate (drops 181 Gerhard/PRJNA512027 + GSE167523;
# adds GSE213621 + GSE162694). Writes NEW artifact corrected_logcpm_c2_nogerhard.rds.
# Canonical corrected_logcpm.rds is NOT touched; consumers NOT repointed.
# ---------------------------------------------------------------------------
set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
Rscript 03c_rebuild_corrected_logcpm_c2_nogerhard.R
