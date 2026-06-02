#!/bin/bash
set -euo pipefail

# =============================================================================
# scRNA-seq Re-integration Pipeline (Phase 2-4)
# Adds GSE136103, re-runs ScaleSC, CellTypist, metadata, pseudobulk DE, concordance
# =============================================================================

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SC_DIR="$BASE/Analysis/SingleCell"
LOGDIR="$SC_DIR/scripts/logs"
mkdir -p "$LOGDIR"

echo "=== scRNA-seq Re-integration Pipeline ==="
echo "Start: $(date)"

# --- Step 1: ScaleSC integration (GPU) ---
JOB1=$(sbatch --parsable \
  --job-name=scalesc_reint \
  --partition=gpu \
  --gres=gpu:1 \
  --cpus-per-task=8 \
  --mem=64G \
  --time=12:00:00 \
  --output="$LOGDIR/scalesc_reint_%j.log" \
  --wrap="bash -c '
set -euo pipefail
eval \"\$(micromamba shell hook --shell=bash)\"
micromamba activate rapids_singlecell
export LD_LIBRARY_PATH=\$(python3 -c \"import nvidia.cusparselt.lib; print(nvidia.cusparselt.lib.__path__[0])\")\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}
cd $BASE
python3 $SC_DIR/integration/scripts/run_scalesc_integration.py --species human --resolution 0.5 --n-hvgs 4000
'")
echo "Step 1: ScaleSC integration -> SLURM $JOB1"

# --- Step 2: CellTypist annotation + pseudobulk export (GPU optional, CPU OK) ---
JOB2=$(sbatch --parsable \
  --dependency=afterok:$JOB1 \
  --job-name=celltypist_reint \
  --partition=cpu \
  --cpus-per-task=8 \
  --mem=96G \
  --time=8:00:00 \
  --output="$LOGDIR/celltypist_reint_%j.log" \
  --wrap="bash -c '
set -euo pipefail
eval \"\$(micromamba shell hook --shell=bash)\"
micromamba activate spatial
cd $BASE
python3 $SC_DIR/integration/scripts/annotate_celltypist.py
'")
echo "Step 2: CellTypist annotation -> SLURM $JOB2 (after $JOB1)"

# --- Step 3: Add sample metadata (condition_harmonized, preparation_method) ---
JOB3=$(sbatch --parsable \
  --dependency=afterok:$JOB2 \
  --job-name=metadata_reint \
  --partition=cpu \
  --cpus-per-task=1 \
  --mem=32G \
  --time=1:00:00 \
  --output="$LOGDIR/metadata_reint_%j.log" \
  --wrap="bash -c '
set -euo pipefail
eval \"\$(micromamba shell hook --shell=bash)\"
micromamba activate spatial
cd $BASE
python3 $SC_DIR/integration/scripts/add_sample_metadata.py
'")
echo "Step 3: Add sample metadata -> SLURM $JOB3 (after $JOB2)"

# --- Step 4: Pseudobulk DE for all cell types ---
JOB4=$(sbatch --parsable \
  --dependency=afterok:$JOB3 \
  --job-name=pbde_reint \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=32G \
  --time=2:00:00 \
  --output="$LOGDIR/pbde_reint_%j.log" \
  --wrap="bash -c '
set -euo pipefail
eval \"\$(micromamba shell hook --shell=bash)\"
micromamba activate rnaseq
cd $BASE
Rscript $SC_DIR/scripts/pseudobulk_de.R
'")
echo "Step 4: Pseudobulk DE -> SLURM $JOB4 (after $JOB3)"

# --- Step 5: Post-GSE136103 concordance (Phase 4.1) ---
JOB5=$(sbatch --parsable \
  --dependency=afterok:$JOB4 \
  --job-name=concordance_post \
  --partition=cpu \
  --cpus-per-task=1 \
  --mem=16G \
  --time=1:00:00 \
  --output="$LOGDIR/concordance_post_%j.log" \
  --wrap="bash -c '
set -euo pipefail
eval \"\$(micromamba shell hook --shell=bash)\"
micromamba activate rnaseq
cd $BASE/RNA-seq
Rscript 36b_pseudobulk_bulk_concordance.R
'")
echo "Step 5: Post-concordance -> SLURM $JOB5 (after $JOB4)"

# --- Step 6: Rebuild atlas with updated sc data (Phase 4.2) ---
JOB6=$(sbatch --parsable \
  --dependency=afterok:$JOB4 \
  --job-name=atlas_rebuild \
  --partition=cpu \
  --cpus-per-task=4 \
  --mem=32G \
  --time=2:00:00 \
  --output="$LOGDIR/atlas_rebuild_%j.log" \
  --wrap="bash -c '
set -euo pipefail
eval \"\$(micromamba shell hook --shell=bash)\"
micromamba activate rnaseq
cd $BASE/RNA-seq
Rscript 45a_integrate_all_sources.R
'")
echo "Step 6: Atlas rebuild -> SLURM $JOB6 (after $JOB4)"

echo ""
echo "=== Pipeline submitted ==="
echo "Step 1 (ScaleSC):    $JOB1"
echo "Step 2 (CellTypist): $JOB2"
echo "Step 3 (Metadata):   $JOB3"
echo "Step 4 (PB DE):      $JOB4"
echo "Step 5 (Concordance):$JOB5"
echo "Step 6 (Atlas):      $JOB6"
echo ""
echo "Monitor: squeue -u $USER -j $JOB1,$JOB2,$JOB3,$JOB4,$JOB5,$JOB6"
