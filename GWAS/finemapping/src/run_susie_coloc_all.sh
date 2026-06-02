#!/bin/bash
# run_susie_coloc_all.sh — Submit SuSiE-COLOC (Script 06) for all 28 GWAS
# Usage: bash src/run_susie_coloc_all.sh
#
# Each GWAS gets a 22-chromosome array job.
# After all complete, run: Rscript src/07_combine_susie_coloc.R

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping

# All 28 GWAS from gwas_registry.tsv
GWAS_LIST=(
  # EUR NAFLD case-control (6)
  2019_31311600_NAFLD_EUR
  2020_32298765_NAFLD_EUR
  2021_34841290_NAFLD_EUR
  2023_36280732_NAFLD_deCode_EUR
  2023_36280732_NAFLD_Intermountain_EUR
  2023_36280732_NAFLD_UKBB_EUR
  # EUR PDFF quantitative (3)
  2021_34128465_PDFF_EUR
  2021_34957434_PDFF_EUR
  2022_36402844_PDFF_EUR
  # EUR liver enzymes (3)
  UKBB_ALT
  UKBB_AST
  UKBB_GGT
  # EUR FinnGen (3)
  FinnGen_NAFLD
  FinnGen_NASH
  FinnGen_HCC
  # EUR cirrhosis/HCC (2)
  Ghouse_Cirrhosis
  Ghouse_HCC
  # EAS case-control (2)
  2020_32514122_Cirrhosis_EAS
  2020_32514122_HCC_EAS
  # EAS liver enzymes (3)
  BBJ_ALT
  BBJ_AST
  BBJ_GGT
  # AFR liver enzymes (3)
  PanUKBB_AFR_ALT
  PanUKBB_AFR_AST
  PanUKBB_AFR_GGT
  # SAS liver enzymes (3)
  PanUKBB_CSA_ALT
  PanUKBB_CSA_AST
  PanUKBB_CSA_GGT
)

echo "Submitting SuSiE-COLOC for ${#GWAS_LIST[@]} GWAS studies"
echo "Each with 22 chromosomes (array jobs)"
echo ""

JOB_IDS=()
for gwas in "${GWAS_LIST[@]}"; do
  JID=$(GWAS_NAME="${gwas}" sbatch --parsable src/06_susie_coloc.sh)
  echo "  ${gwas}: Job ${JID}"
  JOB_IDS+=("${JID}")
done

echo ""
echo "All ${#GWAS_LIST[@]} array jobs submitted (${#JOB_IDS[@]} total)"
echo "Monitor: squeue --me -n susie_coloc"
echo ""
echo "After all complete, run:"
echo "  Rscript src/07_combine_susie_coloc.R"
