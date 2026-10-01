#!/bin/bash
#SBATCH --job-name=agp_ad2_mpra_fit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad2_mpra_fit_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad2_mpra_fit_%j.err

set -uo pipefail

export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export HF_HUB_OFFLINE=1

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
BENCH=$PROJ/Analysis/MASLD_Model_Benchmark
SRC=$PROJ/scripts/analysis/alphagenome_program
OUT=$PROJ/GWAS/finemapping/results/alphagenome_program/ad-arm2-mpra-20260915T114407Z

# The published comparator's own interpreter, so numpy/scipy/sklearn match the arm this contrast is
# nested inside. The reproduction gate inside the script fails the job if it does not return the
# published 0.21622332457761556 and 0.13421551707194201 digit for digit.
PYTHON=$BENCH/executions/environments/transcriptformer-torch2.5.1-21070738/env/bin/python
export PYTHONPATH=$BENCH:$BENCH/src

mkdir -p "$OUT/tables" "$OUT/env"
echo "=== node $(hostname) job ${SLURM_JOB_ID} $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
"$PYTHON" --version 2>&1 | tee "$OUT/env/python_cpu.version.txt"
"$PYTHON" -m pip freeze > "$OUT/env/pip_freeze_cpu.txt" 2>/dev/null

CARD=$(python3 -c "
import json,sys
d=json.load(open('$OUT/embeddings/extraction_receipt_2048bp.json'))
print(d['gpu_card'])
")
echo "card recorded by the extraction: $CARD"

"$PYTHON" "$SRC/i1_fit_reporter_probe.py" \
    --embeddings-2048 "$OUT/embeddings/alphagenome_frozen_trunk_2048bp.npz" \
    --embeddings-16384 "$OUT/embeddings/alphagenome_frozen_trunk_16384bp.npz" \
    --card "$CARD" \
    --out "$OUT"
RC=$?
echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) rc=${RC} ==="
exit $RC
