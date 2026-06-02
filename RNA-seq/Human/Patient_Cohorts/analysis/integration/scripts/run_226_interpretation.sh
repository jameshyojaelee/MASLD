#!/bin/bash
# Submit 226a-d plasma interpretation pipeline with SLURM dependencies
# 226a is the bottleneck (~4h); 226b/c/d run in parallel after it

SCRIPTS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
LOGS="$SCRIPTS/logs"
mkdir -p "$LOGS"

submit_python() {
    local NAME="$1" SCRIPT="$2" CPUS="$3" MEM="$4" TIME="$5" DEP="$6"
    local TMPF=$(mktemp /tmp/sbatch_${NAME}_XXXXXX.sh)
    cat > "$TMPF" << INNEREOF
#!/bin/bash
#SBATCH --job-name=$NAME
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH --output=$LOGS/${NAME}_%j.out
#SBATCH --error=$LOGS/${NAME}_%j.err
cd $SCRIPTS
eval "\$(micromamba shell hook --shell bash)"
micromamba activate spatial
export PYTHONUNBUFFERED=1
python -u $SCRIPT
INNEREOF
    if [ -n "$DEP" ]; then
        sbatch --parsable --dependency=afterok:$DEP "$TMPF"
    else
        sbatch --parsable "$TMPF"
    fi
    rm -f "$TMPF"
}

submit_r() {
    local NAME="$1" SCRIPT="$2" CPUS="$3" MEM="$4" TIME="$5" DEP="$6"
    local TMPF=$(mktemp /tmp/sbatch_${NAME}_XXXXXX.sh)
    cat > "$TMPF" << INNEREOF
#!/bin/bash
#SBATCH --job-name=$NAME
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=$CPUS
#SBATCH --mem=$MEM
#SBATCH --time=$TIME
#SBATCH --output=$LOGS/${NAME}_%j.out
#SBATCH --error=$LOGS/${NAME}_%j.err
cd $SCRIPTS
eval "\$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
Rscript $SCRIPT
INNEREOF
    sbatch --parsable --dependency=afterok:$DEP "$TMPF"
    rm -f "$TMPF"
}

# --- 226a: SHAP + panel size curves (Python, ~4h) ---
JOB_A=$(submit_python "226a_shap" "226a_plasma_shap_panels.py" 16 "64G" "48:00:00" "")
echo "226a submitted: $JOB_A"

# --- 226b: Tissue-plasma bridge (Python, ~10m) — depends on 226a ---
JOB_B=$(submit_python "226b_bridge" "226b_tissue_plasma_bridge.py" 4 "16G" "01:00:00" "$JOB_A")
echo "226b submitted: $JOB_B (after $JOB_A)"

# --- 226c: Cross-platform replication (Python, ~2h) — depends on 226a for protein lists ---
JOB_C=$(submit_python "226c_xplat" "226c_cross_platform_replication.py" 8 "32G" "04:00:00" "$JOB_A")
echo "226c submitted: $JOB_C (after $JOB_A)"

# --- 226d: Etiology pathways (R, ~30m) — depends on 226a for SHAP NPZ ---
JOB_D=$(submit_r "226d_paths" "226d_etiology_pathways.R" 4 "32G" "02:00:00" "$JOB_A")
echo "226d submitted: $JOB_D (after $JOB_A)"

echo ""
echo "Pipeline submitted. 226b/c/d will start automatically when 226a completes."
echo "226e (figures) and 226f (summary) to be submitted after 226a-d finish."
echo "Monitor: squeue -u \$USER | grep 226"
