#!/bin/bash
#SBATCH --job-name=agp_b3-saturation
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_b3-saturation-%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_b3-saturation-%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2

PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
CODE=$PROJ/scripts/analysis/alphagenome_program
OUT=${OUT:?set OUT}
SNAP=$OUT/prespec/region_dir_snapshot_20260914T193225Z.txt

echo "host $(hostname)  start $(date -u +%FT%TZ)  job ${SLURM_JOB_ID:-none}"
sha256sum $OUT/prespec/PRESPECIFICATION.md
echo "--- prespec hash on record ---"; cat $OUT/prespec/PRESPECIFICATION.sha256

echo "=== U2 availability check (spec section 4: no peak-level claim until U2 exists) ==="
SAT=$PROJ/GWAS/finemapping/results/alphagenome_atlas/p4-saturation-20260909T191917Z/raw/saturation
for d in $SAT/*/; do printf "%s\t%s\n" "$(basename $d)" "$(ls -U $d | wc -l)"; done \
  | tee $OUT/tables/universe_directories_on_disk.tsv

echo "=== step 1: fix the sample, then stream region summaries ==="
$PY $CODE/b3_01_sample_and_summarize.py --out $OUT --snapshot $SNAP --workers 8

echo "=== step 2: block bootstrap and associations ==="
$PY $CODE/b3_02_analysis.py --out $OUT

echo "=== provenance ==="
$PY -m pip freeze > $OUT/pip_freeze.txt
$PY -c "import sys,numpy,pandas,h5py; print(sys.version); print('numpy',numpy.__version__,'pandas',pandas.__version__,'h5py',h5py.__version__)" > $OUT/python_versions.txt
{
  printf "path\tsha256\tbytes\trole\n"
  for f in "$PROJ/GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/tables/region_covariates.tsv.gz" \
           "$SAT/../../tables/saturation_universe_summary.json" \
           "$SAT/U1_h3k27ac/chr10_100005812-100006276/request.json" \
           "$CODE/b3_01_sample_and_summarize.py" "$CODE/b3_02_analysis.py" \
           "$OUT/prespec/PRESPECIFICATION.md" "$OUT/prespec/region_dir_snapshot_20260914T193225Z.txt" \
           "$OUT/prespec/sample_regions.tsv"; do
    if [ -f "$f" ]; then printf "%s\t%s\t%s\tinput\n" "$f" "$(sha256sum "$f" | cut -d' ' -f1)" "$(stat -c%s "$f")"; fi
  done
  printf "%s\tNOT_HASHED_TOO_MANY_FILES\t%s\tarchive\n" "$SAT/U1_h3k27ac" "$(ls -U $SAT/U1_h3k27ac | wc -l)"
} > $OUT/MANIFEST.tsv

echo "done $(date -u +%FT%TZ)"
