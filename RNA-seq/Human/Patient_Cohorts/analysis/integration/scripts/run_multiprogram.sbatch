#!/bin/bash
# run_multiprogram.sbatch — Orchestrate multi-program prediction pipeline
# Usage: bash run_multiprogram.sbatch

set -euo pipefail

SCRIPTS=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
LOGDIR=$SCRIPTS/logs
mkdir -p $LOGDIR

echo "=== Multi-Program Prediction Pipeline ==="
echo "Started: $(date)"

# Step 0: NAS subscore recovery (login node, fast)
echo "Running Script 189 (NAS subscore recovery)..."
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
Rscript $SCRIPTS/189_recover_nas_subscores.R

# Step 1: Pseudotime transition signatures
JOB1=$(sbatch --parsable --job-name=191_pseudotime_sig \
  --partition=cpu --cpus-per-task=8 --mem=64G --time=48:00:00 \
  --output=$LOGDIR/191_%j.out --error=$LOGDIR/191_%j.err \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rnaseq && Rscript $SCRIPTS/191_pseudotime_signatures.R'")
echo "Submitted 191: Job $JOB1"

# Step 2: Main multi-program predictor (after pseudotime scores ready)
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 --job-name=190_multiprogram \
  --partition=cpu --cpus-per-task=16 --mem=120G --time=48:00:00 \
  --output=$LOGDIR/190_%j.out --error=$LOGDIR/190_%j.err \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate spatial && cd $SCRIPTS && python 190_multiprogram_predictor.py'")
echo "Submitted 190: Job $JOB2 (depends on $JOB1)"

# Step 3: Genetics-only prediction (independent of pseudotime)
JOB3=$(sbatch --parsable --job-name=192_genetics \
  --partition=cpu --cpus-per-task=8 --mem=64G --time=48:00:00 \
  --output=$LOGDIR/192_%j.out --error=$LOGDIR/192_%j.err \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate spatial && cd $SCRIPTS && python 192_genetics_only_prediction.py'")
echo "Submitted 192: Job $JOB3"

# Step 4: Figures (after 190 and 192 complete)
JOB4=$(sbatch --parsable --dependency=afterok:$JOB2:$JOB3 --job-name=193_figures \
  --partition=cpu --cpus-per-task=4 --mem=32G --time=48:00:00 \
  --output=$LOGDIR/193_%j.out --error=$LOGDIR/193_%j.err \
  --wrap="bash -c 'eval \"\$(micromamba shell hook --shell bash)\" && micromamba activate rnaseq && Rscript $SCRIPTS/193_multiprogram_figures.R'")
echo "Submitted 193: Job $JOB4 (depends on $JOB2, $JOB3)"

echo ""
echo "Pipeline submitted:"
echo "  191 (pseudotime signatures): $JOB1"
echo "  190 (main predictor):        $JOB2 (after $JOB1)"
echo "  192 (genetics-only):         $JOB3 (independent)"
echo "  193 (figures):               $JOB4 (after $JOB2, $JOB3)"
echo ""
echo "Monitor: squeue -u $(whoami) --name=191_pseudotime_sig,190_multiprogram,192_genetics,193_figures"
