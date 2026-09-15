#!/usr/bin/env bash
#SBATCH --job-name=agp_f2-haplotype-v2-report
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f2-haplotype-v2-report-%j.out
set -euo pipefail
export PYTHONNOUSERSITE=1
export F2V2_OUT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping/results/alphagenome_program/f2-haplotype-v2-20260915T063000Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/analysis/alphagenome_program
"$PY" f2v2_03_recompute.py
"$PY" f2v2_04_write_results.py
"$PY" f2v2_07_mark_superseded.py
"$PY" f2v2_05_manifest.py
"$PY" -m pip freeze > "$F2V2_OUT_ROOT/pip_freeze.txt" 2>/dev/null || true
mkdir -p "$F2V2_OUT_ROOT/code_snapshot" "$F2V2_OUT_ROOT/logs"
cp f2_recipe.py f2v2_00_prespec.py f2v2_00_prespec.json f2v2_00_amendment_01.py f2v2_00_amendment_01.json \
   f2v2_01_draw_genespan_null.py f2v2_02_score_genespan_null.py f2v2_02b_loader_equivalence.py \
   f2v2_03_recompute.py f2v2_04_write_results.py f2v2_05_manifest.py f2v2_07_mark_superseded.py \
   f2v2_02_score.sbatch f2v2_02b_probe.sbatch f2v2_06_final.sbatch "$F2V2_OUT_ROOT/code_snapshot/"
cp /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f2-haplotype-v2-*.out "$F2V2_OUT_ROOT/logs/" || true
sacct -u "$USER" -S 2026-09-15 --name=agp_f2-haplotype-v2,agp_f2-haplotype-v2-probe,agp_f2-haplotype-v2-report -P \
  --format=JobID,JobName%30,State,Elapsed,ReqCPUS,ReqMem,Partition,MaxRSS -X > "$F2V2_OUT_ROOT/logs/job_states.tsv" || true
echo F2V2_REPORT_DONE
