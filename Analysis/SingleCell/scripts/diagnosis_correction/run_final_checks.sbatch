#!/bin/bash
#SBATCH --job-name=sc
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/final_checks_%j.out
#SBATCH --error=Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/logs/final_checks_%j.err
set -euo pipefail
TASK="$PWD/Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
test -e "$TASK/MODEL_CHECKS_COMPLETED"
test -e "$TASK/manuscript/REVIEW_STATE.txt"
python3 scripts/manuscript/validate_resource_scope.py > "$TASK/resource_scope_check_final.txt"
python3 - <<'PY' > "$TASK/static_checks_final.txt"
import ast
import csv
import subprocess
from pathlib import Path

root = Path.cwd()
task = root / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
with (task / "modified_files.tsv").open() as handle:
    paths = {Path(row["path"]) for row in csv.DictReader(handle, delimiter="\t")}
paths.update((root / "Analysis/SingleCell/scripts/diagnosis_correction").glob("*"))
paths.update(task / "refit_scripts" / name for name in ["h34.py", "n3_summary.R"])
counts = {"Python": 0, "R": 0, "shell": 0}
rs = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"
for path in sorted(paths):
    if path.suffix == ".py":
        ast.parse(path.read_text(), filename=str(path))
        counts["Python"] += 1
    elif path.suffix == ".R":
        subprocess.run([rs, "-e", "for (p in commandArgs(TRUE)) invisible(parse(p))", str(path)], check=True)
        counts["R"] += 1
    elif path.suffix in {".sbatch", ".sh"}:
        subprocess.run(["bash", "-n", str(path)], check=True)
        counts["shell"] += 1
tracked = subprocess.run(["git", "ls-files"], text=True, capture_output=True, check=True).stdout.splitlines()
relative = [str(path.relative_to(root)) for path in sorted(paths) if str(path.relative_to(root)) in tracked]
subprocess.run(["git", "diff", "--check", "--", *relative], check=True)
print("Syntax and scoped whitespace checks passed:", counts)
PY
touch "$TASK/FINAL_CHECKS_COMPLETED"
