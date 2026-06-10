#!/bin/bash
#SBATCH --job-name=dream
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=%x_lvqw0506_%j.out
#SBATCH --error=%x_lvqw0506_%j.out
# ---------------------------------------------------------------------------
# LVQW staging run — Module 05 (mashr) then 06 (classify), chained after 04.
# Reads + writes lvqw_staging/ ONLY (SEXV3_OUT override). V2_F/CHRX_F canonical.
# ---------------------------------------------------------------------------
# NOTE: NO `set -u` — micromamba/conda activate.d hooks (binutils ADDR2LINE)
# reference unset vars and abort under nounset. Keep -e + pipefail only.
set -eo pipefail
SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/sex_v3
export SEXV3_OUT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/lvqw_staging
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

cd "$SCRIPTS"
# Env split MATCHES the canonical pipeline (the profile.d/micromamba.sh path does
# NOT exist; binary is in ~/.local/bin):
#   - Module 05 (mashr) -> env `motifbreakr` (has mashr; canonical run_module_05.sbatch)
#   - Module 06 (classify) -> env `rnaseq` (data.table + jsonlite; canonical run_module_06.sbatch)
# rnaseq does NOT have mashr, which is why the prior single-env run failed in 4s.
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"

echo "=== LVQW staging Module 05 (mashr; env=motifbreakr) ==="
echo "host: $(hostname)  SEXV3_OUT: $SEXV3_OUT"
echo "started 05: $(date)"
micromamba activate motifbreakr
Rscript 05_mashr_fit.R
micromamba deactivate
echo "finished 05: $(date)"

echo "=== LVQW staging Module 06 (classify; env=rnaseq) ==="
echo "started 06: $(date)"
micromamba activate rnaseq
Rscript 06_classify_v3.R
micromamba deactivate
echo "finished 06: $(date)"
