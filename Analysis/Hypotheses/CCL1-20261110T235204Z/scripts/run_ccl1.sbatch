#!/bin/bash
#SBATCH --job-name=CCL1ident
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CCL1-20261110T235204Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CCL1-20261110T235204Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CCL1-20261110T235204Z/logs/%j.err
set -euo pipefail
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CCL1-20261110T235204Z
OUTD=$RUN/out_ccl1_job${SLURM_JOB_ID}
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python -I "$RUN/scripts/ccl1_identity.py" "$OUTD"
sha256sum "$RUN/scripts/ccl1_identity.py" "$RUN/scripts/run_ccl1.sbatch" "$RUN/PRESPEC_CCL1.md" > "$OUTD/script_sha256.txt"
sha256sum $ROOT/Analysis/Hypotheses/V5-20261021T122723Z/derived/gse298719_qc_cells.h5ad \
  $ROOT/Analysis/Hypotheses/V5-20261021T122723Z/outputs/donor_qc.tsv \
  $ROOT/Analysis/Hypotheses/DIG1-20261024T113951Z/out_snrna_labels/lineage_per_nucleus.tsv.gz \
  $ROOT/Analysis/Hypotheses/V3-20261008T015239Z/pseudobulk/donor_metadata.tsv \
  $ROOT/Analysis/Hypotheses/V3-20261008T015239Z/downloads_untrusted/GSM*_CB_raw_feature_bc_matrix_filtered.h5 \
  $ROOT/data/GSE202379/metadata/donor_pairing.csv $ROOT/data/GSE244832/metadata/donor_pairing.csv \
  $ROOT/Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/metadata/donor_metadata_stage_corrected_donor.tsv \
  $ROOT/data/public_r4/GSE289173/derived_H90-20261015T120337Z/gse289173_qc_nuclei.h5ad \
  $ROOT/Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad > "$OUTD/input_sha256.txt"
echo done
