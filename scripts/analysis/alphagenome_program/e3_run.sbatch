#!/bin/bash
#SBATCH --job-name=agp_e3-reporter-context
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e3_run_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_e3_run_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
CODE="$PROJ/scripts/analysis/alphagenome_program"
OUT="$1"
PY=/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_atlas/bin/python

echo "host $(hostname) start $(date -u +%Y-%m-%dT%H:%M:%SZ) job ${SLURM_JOB_ID}"

"$PY" "$CODE/e3_score_api.py" "$OUT" --max-calls 1200
echo "scoring done $(date -u +%Y-%m-%dT%H:%M:%SZ)"

"$PY" "$CODE/e3_analyse.py" "$OUT"
echo "analysis done $(date -u +%Y-%m-%dT%H:%M:%SZ)"

"$PY" -m pip freeze > "$OUT/pip_freeze_alphagenome_atlas.txt" 2>/dev/null || true
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze \
    > "$OUT/pip_freeze_rnaseq.txt" 2>/dev/null || true

{
  printf 'path\trole\tsha256\n'
  for f in \
    "$PROJ/Analysis/MASLD_Model_Benchmark/executions/gse281364-validated-reconstruction-21066470/validated/elements.tsv:mpra_elements" \
    "$PROJ/Analysis/MASLD_Model_Benchmark/executions/gse281364-validated-reconstruction-21066470/validated/replicate_outcomes.tsv.gz:mpra_replicate_outcomes" \
    "$PROJ/Analysis/MASLD_Model_Benchmark/executions/gse281364-borzoi-native-fixture-21083008/fixture/manifest.tsv:block_rule_verification" \
    "$PROJ/Analysis/MASLD_Model_Benchmark/executions/gse281364-borzoi-native-fixture-21083008/fixture/long_range_components.tsv:block_rule_verification" \
    "$PROJ/scripts/analysis/alphagenome_program/e3_prepare.py:code" \
    "$PROJ/scripts/analysis/alphagenome_program/e3_score_api.py:code" \
    "$PROJ/scripts/analysis/alphagenome_program/e3_analyse.py:code" \
    "$PROJ/scripts/analysis/alphagenome_atlas/atlas_query.py:quota_retry_wrapper_read_only" \
    "$PROJ/scripts/analysis/alphagenome_atlas/63_indel_rescue.py:api_call_pattern_read_only" \
    "$PROJ/scripts/analysis/alphagenome_program/IMPLEMENTATION_SPEC.md:spec" \
    "$OUT/prespec/PRESPECIFICATION.md:prespecification" ; do
    p="${f%%:*}"; r="${f##*:}"
    printf '%s\t%s\t%s\n' "$p" "$r" "$(sha256sum "$p" | cut -d' ' -f1)"
  done
  printf '%s\t%s\t%s\n' \
    "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa" \
    "reference_fasta_no_sha_too_large" "not_computed"
} > "$OUT/MANIFEST.tsv"

echo "end $(date -u +%Y-%m-%dT%H:%M:%SZ)"
