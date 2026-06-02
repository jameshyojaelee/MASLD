#!/bin/bash
# Submit 225b/225d/225e models as parallel SLURM jobs — one per model (or small group)
#
# Uses a template sbatch file to avoid --wrap shell compatibility issues.

SCRIPTS="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"
cd "$SCRIPTS"

COMMON_SBATCH="--partition=cpu --cpus-per-task=16 --mem=64G --time=48:00:00"

# Template sbatch script generator
submit_model() {
    local SCRIPT="$1"
    local MODEL="$2"
    local SUFFIX="$3"
    local JOBNAME="$4"

    local TMPSCRIPT=$(mktemp /tmp/sbatch_${JOBNAME}_XXXXXX.sh)
    cat > "$TMPSCRIPT" << INNEREOF
#!/bin/bash
#SBATCH --job-name=$JOBNAME
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=$SCRIPTS/logs/${JOBNAME}_%j.out
#SBATCH --error=$SCRIPTS/logs/${JOBNAME}_%j.err
cd $SCRIPTS
eval "\$(micromamba shell hook --shell bash)"
micromamba activate spatial
export MODEL_SUBSET="$MODEL"
export OUT_SUFFIX="$SUFFIX"
python $SCRIPT
INNEREOF
    sbatch "$TMPSCRIPT"
    rm -f "$TMPSCRIPT"
}

echo "=== 225b: remaining models ==="
submit_model "225b_plasma_sweep_boosting.py" "lightgbm" "lgb" "225b_lgb"
submit_model "225b_plasma_sweep_boosting.py" "adaboost" "ada" "225b_ada"

echo ""
echo "=== 225d: remaining ordinal models (11 done, ~23 remaining) ==="
submit_model "225d_plasma_sweep_ordinal.py" "lightgbm" "lgb" "225d_lgb"
submit_model "225d_plasma_sweep_ordinal.py" "linear_svc,rbf_svm,nu_svc" "svm" "225d_svm"
submit_model "225d_plasma_sweep_ordinal.py" "mlp" "mlp" "225d_mlp"
submit_model "225d_plasma_sweep_ordinal.py" "adaboost" "ada" "225d_ada"
submit_model "225d_plasma_sweep_ordinal.py" "catboost" "cat" "225d_cat"
submit_model "225d_plasma_sweep_ordinal.py" "easy_ensemble,balanced_bagging" "imb" "225d_imb"
submit_model "225d_plasma_sweep_ordinal.py" "ordinal_forest,coral" "mord" "225d_mord"
submit_model "225d_plasma_sweep_ordinal.py" "platt_svm,calibrated_ensemble" "cal" "225d_cal"
submit_model "225d_plasma_sweep_ordinal.py" "gaussian_process" "gp" "225d_gp"
submit_model "225d_plasma_sweep_ordinal.py" "kernel_ridge,elasticnet_reg" "reg" "225d_reg"
submit_model "225d_plasma_sweep_ordinal.py" "smote_rf" "smote" "225d_smote"
submit_model "225d_plasma_sweep_ordinal.py" "tabnet,tabpfn" "tab" "225d_tab"
submit_model "225d_plasma_sweep_ordinal.py" "flaml_30s,flaml_60s" "flaml" "225d_flaml"
submit_model "225d_plasma_sweep_ordinal.py" "voting_lr_rf_hgb,stacking_lr_rf_hgb" "ens" "225d_ens"

echo ""
echo "=== 225e: remaining etiology models (3 done, ~29 remaining) ==="
submit_model "225e_plasma_sweep_etiology.py" "sgd_elasticnet,bayesian_ridge" "sgd_bayes" "225e_sgd"
submit_model "225e_plasma_sweep_etiology.py" "random_forest,balanced_rf,extra_trees" "tree" "225e_tree"
submit_model "225e_plasma_sweep_etiology.py" "gradient_boosting" "gbm" "225e_gbm"
submit_model "225e_plasma_sweep_etiology.py" "hist_gradient_boosting" "hgbm" "225e_hgbm"
submit_model "225e_plasma_sweep_etiology.py" "xgboost" "xgb" "225e_xgb"
submit_model "225e_plasma_sweep_etiology.py" "lightgbm" "lgb" "225e_lgb"
submit_model "225e_plasma_sweep_etiology.py" "linear_svc,rbf_svm,nu_svc" "svm" "225e_svm"
submit_model "225e_plasma_sweep_etiology.py" "mlp" "mlp" "225e_mlp"
submit_model "225e_plasma_sweep_etiology.py" "adaboost" "ada" "225e_ada"
submit_model "225e_plasma_sweep_etiology.py" "catboost" "cat" "225e_cat"
submit_model "225e_plasma_sweep_etiology.py" "easy_ensemble,balanced_bagging" "imb" "225e_imb"
submit_model "225e_plasma_sweep_etiology.py" "platt_svm,calibrated_ensemble" "cal" "225e_cal"
submit_model "225e_plasma_sweep_etiology.py" "gaussian_process" "gp" "225e_gp"
submit_model "225e_plasma_sweep_etiology.py" "kernel_ridge,elasticnet_reg" "reg" "225e_reg"
submit_model "225e_plasma_sweep_etiology.py" "smote_rf" "smote" "225e_smote"
submit_model "225e_plasma_sweep_etiology.py" "tabnet,tabpfn" "tab" "225e_tab"
submit_model "225e_plasma_sweep_etiology.py" "flaml_30s,flaml_60s" "flaml" "225e_flaml"
submit_model "225e_plasma_sweep_etiology.py" "voting_lr_rf_hgb,stacking_lr_rf_hgb" "ens" "225e_ens"

echo ""
echo "Total jobs submitted. Use 'squeue -u \$USER | grep 225' to monitor."
