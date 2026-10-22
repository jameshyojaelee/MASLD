#!/bin/bash
#SBATCH --job-name=7601
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F6/logs/7601_%j.log
set -euo pipefail
W=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F6
P=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/review-r3-fixes-20261015T112148Z/F6/repair
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 TMPDIR=$W/tmp MPLCONFIGDIR=$W/tmp
PYS=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
R=$W/repair_20261015T172841Z
TS=$(date -u +%Y%m%dT%H%M%SZ)
O=$P/panels_$TS
mkdir $O
$PYS $W/f6_render_panels_repair.py $R $O
cp $R/R15_03_gene_effect_sidecar_repair.tsv.gz $R/R15_05_assay_applicability_corrected_repair.tsv.gz $R/A_R07_02_seven_susie_untestable_genes.tsv $R/B_R07_02_pph3_tier12.tsv $R/C_R07_05_bulk_ids_absent.tsv $R/C_R07_05_t12_ids_absent.tsv $R/F_R10_03_myojin_frozen_code_refit.tsv $R/F_R10_03_permutation_audit.tsv $R/f6_repair_summary.json $O/
cp $W/f6_repair.py $W/f6_render_panels_repair.py $W/phase_c_stats_copy.py $O/
$PYS -m pip freeze > $O/pip_freeze_spatial.txt
cp $R/pip_freeze_rnaseq.txt $O/
for p in $O/*.pdf; do pdftoppm -r 200 -png -singlefile $p ${p%.pdf}; pdftotext -layout $p ${p%.pdf}.txt; done
echo DONE $O
