#!/bin/bash
#SBATCH --job-name=N3source_verify
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/N3-source-diagnosis-20261014/logs/%x-%j.out
set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
TASK_ROOT=$PWD/Analysis/Hypotheses/N3-source-diagnosis-20261014
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 - <<'PY'
import csv,hashlib
from pathlib import Path
b=Path('Analysis/Hypotheses/N3-source-diagnosis-20261014')
with (b/'per_library_input_sha256.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['input','sha256'])
 for p in sorted((b/'corrected/per_library').glob('*.rds')):
  w.writerow([str(p),hashlib.sha256(p.read_bytes()).hexdigest()])
PY
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript "$TASK_ROOT/scripts/04_ambient_summary.R" "$TASK_ROOT/corrected" "$TASK_ROOT/raw"
python3 "$TASK_ROOT/scripts/verify_terminal.py"
