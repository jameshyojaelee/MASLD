#!/bin/bash
#SBATCH --job-name=conc_lfsr
#SBATCH --partition=io
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/patient_concordance_lfsr_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/patient_concordance_lfsr_%j.err

# lfsr variant of the patient-concordance figure family.
# Same analytical body as the raw-LFC scripts; DEG set = ashr lfsr<0.05 (sign-confidence),
# effect = shrunk_logFC. Output -> figures/supplementary/figS_methods_validation/lfc_sensitivity/patient_concordance_lfsr/
set -eo pipefail
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

mkdir -p figures/supplementary/figS_methods_validation/lfc_sensitivity/patient_concordance_lfsr

echo "=== [1/2] figS_patient_concordance_lfsr.R (A-D sweep + density/robustness panels) ==="
Rscript scripts/figures/figS_patient_concordance_lfsr.R
echo "=== [2/2] figS_pi_concordance_panels_lfsr.R (E1-E4 stage/waterfall/ROC/cohort) ==="
Rscript scripts/figures/figS_pi_concordance_panels_lfsr.R
echo "=== DONE ==="
ls -1 figures/supplementary/figS_methods_validation/lfc_sensitivity/patient_concordance_lfsr/*.pdf | wc -l
