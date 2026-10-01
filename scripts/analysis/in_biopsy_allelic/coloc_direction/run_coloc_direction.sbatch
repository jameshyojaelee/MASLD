#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=90:00:00
# B-COLOC direction table: self-tests, then coloc_direction.py.
# Usage: [COMPARE_TO=<previous run dir>] sbatch --job-name=<n> --output=<exec>/logs/%j_coloc_direction.out \
#          run_coloc_direction.sbatch <out_dir> [--studies A,B --chroms 19,22]
# Reads no MASLD allelic data; outputs are per-gene / per-variant summaries only.
set -euo pipefail
out=$1; shift
here=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/in_biopsy_allelic/coloc_direction
py="env PYTHONNOUSERSITE=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python"
module load bcftools/1.21 htslib/1.21 >/dev/null 2>&1
echo "host $(hostname) job ${SLURM_JOB_ID:-none} start $(date -u +%FT%TZ)"
bcftools --version | head -1
$py -m pytest -q -p no:cacheprovider "$here/tests/test_direction_rule.py"
work=$out/work   # lead-position extracts of public GWAS/eQTL/1kGP files, kept for re-derivation
$py "$here/coloc_direction.py" --out-dir "$out" --work-dir "$work" \
  --liftover-bin "$here/liftover.sh" --threads "${SLURM_CPUS_PER_TASK:-8}" "$@"
if [ -n "${COMPARE_TO:-}" ]; then
  $py "$here/compare_runs.py" "$COMPARE_TO" "$out" "$out/changes_vs_previous.tsv"
fi
$py -m pip freeze > "$out/pip_freeze.txt"
sha256sum "$here/coloc_direction.py" "$here/eqtl_conditional_sign.R" "$here/tests/test_direction_rule.py" "$here/run_coloc_direction.sbatch" "$here/liftover.sh" "$here/compare_runs.py" > "$out/code_sha256.txt"
echo "end $(date -u +%FT%TZ)"
