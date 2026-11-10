#!/usr/bin/env bash
#SBATCH --job-name=masld_atac_da
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/full_da_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods/logs/full_da_%j.log

set -euo pipefail
TASK_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/ATAC/Regulatory_Methods
INPUT=${1:?Usage: sbatch run_full_da.sbatch INPUT_ROOT NEW_OUTPUT}
OUTPUT=${2:?A new output directory is required}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=1
for COHORT in GSE244832 GSE281367; do
  for LINEAGE in hepatocyte stellate macrophage; do
    test -f "$INPUT/count/$COHORT/$LINEAGE/diagnostics.json"
    test -f "$INPUT/binary_counts/$COHORT/$LINEAGE.mtx.gz"
    test -f "$INPUT/balanced_counts/$COHORT/$LINEAGE.mtx.gz"
    test -f "$INPUT/balanced_binary_counts/$COHORT/$LINEAGE.mtx.gz"
  done
done
python3 - "$INPUT" <<'PY'
import json
import sys
from pathlib import Path
root = Path(sys.argv[1])
for cohort in ("GSE244832", "GSE281367"):
    for lineage in ("hepatocyte", "stellate", "macrophage"):
        path = root / "count" / cohort / lineage / "diagnostics.json"
        result = json.loads(path.read_text())
        if result.get("pilot") is not False or result.get("v3_conservation") != "exact_all_donor_peak_cells":
            raise ValueError(f"Full counts have not passed exact conservation: {path}")
        if result.get("cohort") != cohort or result.get("lineage") != lineage:
            raise ValueError(f"Count identity mismatch: {path}")
PY
test ! -e "$OUTPUT"
mkdir -p "$OUTPUT"
cp "$TASK_ROOT/config/analysis.json" "$OUTPUT/analysis.json"
cp "$TASK_ROOT/scripts/compare_da.R" "$OUTPUT/compare_da.R"
cp "$TASK_ROOT/scripts/render_results.R" "$OUTPUT/render_results.R"
cp "$TASK_ROOT/scripts/compare_da_sensitivity.R" "$OUTPUT/compare_da_sensitivity.R"
sha256sum "$OUTPUT/analysis.json" "$OUTPUT/compare_da.R" "$OUTPUT/render_results.R" "$OUTPUT/compare_da_sensitivity.R" > "$OUTPUT/source_sha256.txt"
RNA_R=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
/usr/bin/time -v "$RNA_R" "$OUTPUT/compare_da.R" --config "$OUTPUT/analysis.json" \
  --output "$OUTPUT/full" --binary-root "$INPUT/binary_counts" \
  --workers "$SLURM_CPUS_PER_TASK" --null-reps 1000
/usr/bin/time -v "$RNA_R" "$OUTPUT/compare_da.R" --config "$OUTPUT/analysis.json" \
  --output "$OUTPUT/cell_downsample" --counts-root "$INPUT/balanced_counts" \
  --binary-root "$INPUT/balanced_binary_counts" --no-baseline-check \
  --sensitivity-label cell_downsample --workers "$SLURM_CPUS_PER_TASK" --null-reps 1000
"$RNA_R" "$OUTPUT/render_results.R" --da "$OUTPUT/full" --output "$OUTPUT/figures_full"
"$RNA_R" "$OUTPUT/render_results.R" --da "$OUTPUT/cell_downsample" --output "$OUTPUT/figures_cell_downsample"
"$RNA_R" "$OUTPUT/compare_da_sensitivity.R" --full "$OUTPUT/full" \
  --downsample "$OUTPUT/cell_downsample" --output "$OUTPUT/common_peak_cell_sensitivity"
printf 'All accessibility comparisons and PDF panels completed. Candidate only.\n' > "$OUTPUT/COMPLETE.txt"
