#!/bin/bash
#SBATCH --job-name=A9_chrx_aneup
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/A9_chrx_aneup_%j.out
#SBATCH --error=logs/A9_chrx_aneup_%j.err

# Team A, agent A9: chrX/XCI stratification (T2) + aneuploidy detection (T8)
# Date: 2026-05-11

set -eo pipefail   # NOT -u (micromamba activation uses unbound vars)

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
INTEG=$PROJ/RNA-seq/Human/Patient_Cohorts/analysis/integration
SCRIPTS=$INTEG/scripts

mkdir -p $SCRIPTS/logs

cd $SCRIPTS

# Activate rnaseq env (R + rtracklayer + edgeR + data.table + ggplot2)
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

echo "== Env =="
echo "Python: $(which python)"
echo "Rscript: $(which Rscript)"
echo "Rscript version: $(Rscript --version 2>&1 | head -1)"
echo "Working dir: $(pwd)"

START=$(date +%s)

echo
echo "===================================================="
echo "==  T2: 14.4g_chrx_stratification_T2.R          =="
echo "===================================================="
Rscript --vanilla 14.4g_chrx_stratification_T2.R 2>&1
T2_RC=${PIPESTATUS[0]}
echo "T2 return code: $T2_RC"

echo
echo "===================================================="
echo "==  T8: 14.4h_aneuploidy_T8.R                   =="
echo "===================================================="
Rscript --vanilla 14.4h_aneuploidy_T8.R 2>&1
T8_RC=${PIPESTATUS[0]}
echo "T8 return code: $T8_RC"

END=$(date +%s)
echo
echo "Elapsed: $((END - START)) sec"
echo "Outputs:"
ls -la $PROJ/RNA-seq/results/audit_sensitivity/sex_xci/ 2>&1 || echo "(sex_xci missing)"
ls -la $PROJ/RNA-seq/results/audit_sensitivity/sex_aneuploidy/ 2>&1 || echo "(sex_aneuploidy missing)"
ls -la $PROJ/figures/supplementary/figS_xci_stratification.pdf 2>&1 || echo "(figS_xci missing)"

exit $(( T2_RC + T8_RC ))
