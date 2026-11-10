#!/usr/bin/env bash
# One original outer fold per job; donor assignments and scientific settings persist.
#SBATCH --job-name=masld_atac_fold
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/network_fold_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/network_fold_%j.log
set -euo pipefail
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods
SOURCE_ROOT=${REGULATORY_SOURCE_ROOT:-"$TASK_ROOT"}
INPUT=${1:?Usage: sbatch run_network_fold.sbatch INPUT GENES NEW_OUTPUT CONFIG [RNA_INPUT]}
GENES=${2:?Gene annotation required}
OUTPUT=${3:?New output required}
CONFIG=${4:?Explicit single-fold configuration required}
RNA_INPUT=${5:-}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 - "$INPUT" "$RNA_INPUT" "$CONFIG" <<'PY'
import json
import sys
from pathlib import Path
root = Path(sys.argv[1])
if (root / "diagnostics.json").exists():
    result = json.loads((root / "diagnostics.json").read_text())
    if result.get("pilot") is not False or result.get("v3_conservation") != "exact_all_donor_peak_cells":
        raise ValueError("Disease input must pass full count conservation")
elif (root / "paired_export.json").exists():
    result = json.loads((root / "paired_export.json").read_text())
    if result.get("mode") != "full" or not sys.argv[2] or not (Path(sys.argv[2]) / "RNA_export.json").exists():
        raise ValueError("Paired input requires qualified full raw RNA export")
else:
    raise ValueError("Missing qualified input receipt")
config = json.loads(Path(sys.argv[3]).read_text())
folds = config["network"].get("selected_outer_folds", [])
if len(folds) != 1 or folds[0] not in range(5):
    raise ValueError("Exactly one original outer fold 0–4 is required")
PY
test ! -e "$OUTPUT"
mkdir -p "$OUTPUT/source"
cp "$CONFIG" "$OUTPUT/source/analysis.json"
for SCRIPT in network_compare.py network_cicero.R network_directnet_native.R; do
  cp "$SOURCE_ROOT/scripts/$SCRIPT" "$OUTPUT/source/$SCRIPT"
done
sha256sum "$OUTPUT/source/"* > "$OUTPUT/source_sha256.txt"
MODEL_PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
RNA_ARGS=()
if [[ -n "$RNA_INPUT" ]]; then RNA_ARGS=(--rna-root "$RNA_INPUT"); fi
/usr/bin/time -v "$MODEL_PY" "$OUTPUT/source/network_compare.py" \
  --config "$OUTPUT/source/analysis.json" --input-root "$INPUT" --genes "$GENES" \
  --output "$OUTPUT/analysis" --max-genes 0 "${RNA_ARGS[@]}"
python3 - "$OUTPUT/analysis" "$CONFIG" <<'PY'
import json
import sys
from pathlib import Path
root = Path(sys.argv[1])
result = json.loads((root / "shard_completion.json").read_text())
expected = json.loads(Path(sys.argv[2]).read_text())["network"]["selected_outer_folds"]
if result.get("status") != "completed" or result.get("selected_outer_folds") != expected:
    raise ValueError("Outer-fold shard failed qualification")
if json.loads((root / "completion.json").read_text()).get("status") == "completed":
    raise ValueError("A single-fold shard must not claim full five-fold inference")
PY
printf 'Outer fold complete; five-fold merge and scientific review required.\n' > "$OUTPUT/COMPLETE.txt"
