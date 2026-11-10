#!/usr/bin/env bash
# One cohort and lineage per job; submit individually, without an array.
#SBATCH --job-name=masld_atac_links
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/network_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/network_%j.log

set -euo pipefail
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods
SOURCE_ROOT=${REGULATORY_SOURCE_ROOT:-"$TASK_ROOT"}
INPUT=${1:?Usage: sbatch run_network.sbatch INPUT_ROOT GENES_TSV NEW_OUTPUT [RNA_ROOT]}
GENES=${2:?Gene annotation is required}
OUTPUT=${3:?A new output directory is required}
RNA_INPUT=${4:-}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 - "$INPUT" "$RNA_INPUT" <<'PY'
import json
import sys
from pathlib import Path
root = Path(sys.argv[1])
if (root / "diagnostics.json").exists():
    result = json.loads((root / "diagnostics.json").read_text())
    if result.get("pilot") is not False or result.get("v3_conservation") != "exact_all_donor_peak_cells":
        raise ValueError("Disease-cohort input must pass full count conservation")
elif (root / "paired_export.json").exists():
    result = json.loads((root / "paired_export.json").read_text())
    if result.get("mode") != "full" or not sys.argv[2]:
        raise ValueError("Paired-cohort full analysis requires qualified full RNA input")
    if not (Path(sys.argv[2]) / "RNA_export.json").exists():
        raise ValueError("Missing exact raw-RNA export receipt")
else:
    raise ValueError("Missing qualified input receipt")
PY
test ! -e "$OUTPUT"
mkdir -p "$OUTPUT/source"
cp "$SOURCE_ROOT/config/analysis.json" "$OUTPUT/source/analysis.json"
cp "$SOURCE_ROOT/scripts/network_compare.py" "$OUTPUT/source/network_compare.py"
cp "$SOURCE_ROOT/scripts/network_cicero.R" "$OUTPUT/source/network_cicero.R"
cp "$SOURCE_ROOT/scripts/network_directnet_native.R" "$OUTPUT/source/network_directnet_native.R"
cp "$SOURCE_ROOT/scripts/render_results.R" "$OUTPUT/source/render_results.R"
sha256sum "$OUTPUT/source/"* > "$OUTPUT/source_sha256.txt"
MODEL_PY=/gpfs/commons/home/jameslee/micromamba/envs/rapids_singlecell/bin/python
RNA_ARGS=()
if [[ -n "$RNA_INPUT" ]]; then RNA_ARGS=(--rna-root "$RNA_INPUT"); fi
/usr/bin/time -v "$MODEL_PY" "$OUTPUT/source/network_compare.py" \
  --config "$OUTPUT/source/analysis.json" --input-root "$INPUT" --genes "$GENES" \
  --output "$OUTPUT/analysis" --max-genes 0 "${RNA_ARGS[@]}"
python3 - "$OUTPUT/analysis" <<'PY'
import csv
import json
import sys
from pathlib import Path
root = Path(sys.argv[1])
result = json.loads((root / "completion.json").read_text())
if result.get("status") != "completed":
    raise ValueError("Network analysis did not complete")
tables = list(root.glob("*/method_execution_summary.tsv"))
if not tables:
    raise ValueError("No evaluable cohort-lineage network analysis")
for path in tables:
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if int(row["n_failed_stages"]) > 0:
                raise RuntimeError(f"Incomplete method execution: {path}: {row['method']}")
PY
RNA_R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
"$RNA_R" "$OUTPUT/source/render_results.R" --network "$OUTPUT/analysis" --output "$OUTPUT/figures"
printf 'Candidate network comparisons completed; scientific review required.\n' > "$OUTPUT/COMPLETE.txt"
