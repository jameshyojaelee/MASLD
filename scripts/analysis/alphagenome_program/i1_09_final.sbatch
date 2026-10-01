#!/bin/bash
#SBATCH --job-name=agp_ad1_final
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad1_final_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_ad1_final_%j.err

# Zero-shot endpoint, job table, manifest and RESULTS.md. No GPU: the scores exist.
# The zero-shot step is run only when the gate passed; the gate verdict is read from disk here rather
# than passed in, so this script cannot report a zero-shot number the gate did not license.

set -uo pipefail

export PYTHONNOUSERSITE=1
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba

MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
ENV=alphagenome_local
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
SRC=$PROJ/scripts/analysis/alphagenome_program
RUN=$PROJ/GWAS/finemapping/results/alphagenome_program/ad1-gate-20260915T114615Z
LENGTH=${AD1_LENGTH:?set AD1_LENGTH}

echo "=== node $(hostname) job $SLURM_JOB_ID $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="

echo "=== pilot analysis, rerun with the matched-length table ==="
$MM run -n $ENV python "$SRC/i1_03_pilot_analysis.py" --run "$RUN"
echo "  exit=$?"

VERDICT=$($MM run -n $ENV python -c \
  "import json;print(json.load(open('$RUN/tables/ad1_gate.json'))['gate_verdict'])")
echo "gate verdict read from disk: $VERDICT"

if [ "$VERDICT" = "PASS" ]; then
  echo "=== zero-shot endpoint ==="
  $MM run -n $ENV python "$SRC/i1_05_zeroshot.py" \
      --run "$RUN" --scores "$RUN/raw/gate_${LENGTH}.tsv"
  echo "  exit=$?"
else
  echo "gate did not pass; the arm stops here and no zero-shot number is computed"
fi

echo "=== job table ==="
cp -f "$SRC"/i1_*.sbatch "$RUN/sbatch/" 2>/dev/null
cp -f "$SRC"/i1_*.py "$RUN/sbatch/" 2>/dev/null
( cd "$RUN/sbatch" && sha256sum i1_* > analysis_scripts.sha256 )
sacct -j "${AD1_JOBS:?set AD1_JOBS as a comma-separated job id list}" \
  --format=JobID,JobName%24,State,Elapsed,MaxRSS,NodeList%18,AllocTRES%90 -X -P \
  | tr '|' '\t' > "$RUN/tables/jobs.tsv"
cat "$RUN/tables/jobs.tsv"

echo "=== manifest ==="
$MM run -n $ENV python "$SRC/i1_07_manifest.py" --run "$RUN"
echo "  exit=$?"

echo "=== RESULTS.md ==="
$MM run -n $ENV python "$SRC/i1_08_write_results.py" --run "$RUN"
echo "  exit=$?"

echo "=== scope validator ==="
$MM run -n $ENV python "$PROJ/scripts/manuscript/validate_resource_scope.py" || true

echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
