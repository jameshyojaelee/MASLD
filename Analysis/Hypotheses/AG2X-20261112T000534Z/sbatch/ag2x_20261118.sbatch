#!/usr/bin/env bash
#SBATCH --job-name=ag2x_20261118
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/ag2x_20261118_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/AG2X-20261112T000534Z/logs/ag2x_20261118_%j.err
set -euo pipefail

REPO="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
CODE="${REPO}/Analysis/Hypotheses/AG2X-20261112T000534Z/code"
PY="/gpfs/commons/home/jameslee/micromamba/envs/alphagenome2/bin/python"
RS="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
OUT="${REPO}/Analysis/Hypotheses/AG2X-20261112T000534Z/out_20261118"
mkdir -p "${OUT}/stageA" "${OUT}/stageB" "${OUT}/stageB2" "${OUT}/stageD/tables" "${OUT}/stageD/panels"

export PYTHONHASHSEED=20261111
cd "${CODE}"

echo "=== Stage A (20261118) ==="
"$PY" -u "${CODE}/a20_analyze_20261118.py"

echo "=== Stage B (20261118) ==="
"$PY" -u "${CODE}/b20_analyze_stageB_20261118.py"
"$RS" "${CODE}/b21_b3_mh_20261118.R" "${OUT}/stageB"
"$PY" -u "${CODE}/b23_b3a_corrected_20261118.py"
"$RS" "${CODE}/b24_b3a_corrected_mh_20261118.R" "${OUT}/stageB"

echo "=== Stage B2 (20261118) ==="
"$PY" -u "${CODE}/b2_20_analyze_20261118.py"
"$RS" "${CODE}/b2_21_b3_mh_20261118.R" "${OUT}/stageB2"

echo "=== Stage D (20261118) ==="
"$PY" -u "${CODE}/d20_panels_stageD_20261118.py"
"$PY" -u "${CODE}/d21_captions_stageD_20261118.py"
"$PY" -u "${CODE}/d22_checks_stageD_20261118.py"

echo "=== DONE ==="
