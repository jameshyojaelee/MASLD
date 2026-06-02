#!/usr/bin/env bash
# run_all.sh — chained sbatch for mega_validation pipeline.
#
# Parallelism (8 parallel job slots max — most run concurrently):
#   - perstudy refresh + arm1 (edgeR-QL) + arm2 (voomLmFit): kick off in parallel
#   - arm3 (metafor), arm4 (mashr): start after perstudy
#   - concordance: starts after all 4 arms (arm4 with afterany so a mashr
#                  failure does not block the rest)
#   - figS + atlas refresh: parallel after concordance
set -euo pipefail

PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
INT="$PROJECT_ROOT/RNA-seq/Human/Patient_Cohorts/analysis/integration"
SCRIPTS="$INT/scripts/mega_validation"
LOG="$SCRIPTS/logs"
mkdir -p "$LOG"

SKIP_PERSTUDY="${SKIP_PERSTUDY:-0}"
QOS="${QOS:-nslab}"

# Validate prerequisites
for f in \
  "$INT/results/integration/merged_dge.rds" \
  "$INT/results/integration/meta_matched.rds" \
  "$INT/results/integration/dream_results.csv" \
  "$PROJECT_ROOT/config/human_datasets.yaml"
do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: missing prerequisite $f" >&2; exit 1
  fi
done

submit() {
  local name=$1; shift
  local deps=$1; shift           # raw --dependency arg (may be empty)
  local partition=$1; shift
  local cpus=$1; shift
  local mem=$1; shift
  local time=$1; shift
  local script=$1; shift
  local env_name="${1:-rnaseq}"  # optional 8th arg: env name (default rnaseq)
  local dep_flag=""
  if [[ -n "$deps" ]]; then dep_flag="--dependency=$deps"; fi
  # NB: --wrap defaults to /bin/sh on this cluster; force bash via `bash -lc`
  # so micromamba's shell hook works (source is not a POSIX-sh builtin).
  sbatch --parsable \
         --partition=$partition --qos="$QOS" \
         --cpus-per-task=$cpus --mem=$mem --time=$time \
         --job-name="mvL_$name" \
         --output="$LOG/${name}_%j.out" \
         --error="$LOG/${name}_%j.err" \
         $dep_flag \
         --wrap "bash -lc 'set -e; source ~/.bashrc; micromamba activate $env_name; cd \"$PROJECT_ROOT\"; export MASLD_PROJECT_ROOT=\"$PROJECT_ROOT\"; Rscript \"$script\"'"
}

# Partition routing — cpu is currently congested; use io for small jobs and
# bigmem (>=512G per cluster policy) for the two heavy stats arms.
PART_SMALL="${PART_SMALL:-io}"
PART_HEAVY="${PART_HEAVY:-bigmem}"
HEAVY_MEM="${HEAVY_MEM:-512G}"

# --- Tier 1: independent kickoff (parallel) ---
deps_step1=""
jid_perstudy=""
if [[ "$SKIP_PERSTUDY" == "0" ]]; then
  jid_perstudy=$(submit perstudy "" "$PART_SMALL" 8 32G 24:00:00 "$INT/scripts/02_per_study_de.R")
  echo "perstudy:    $jid_perstudy ($PART_SMALL)"
  deps_step1="afterok:$jid_perstudy"
fi

# Arms 1 + 2 do NOT need per_study refresh — read merged_dge.rds directly.
jid_eql=$(submit eql "" "$PART_SMALL" 8 64G 48:00:00 "$SCRIPTS/01_edgeRql_cohortFE.R")
jid_vlm=$(submit vlm "" "$PART_SMALL" 8 64G 48:00:00 "$SCRIPTS/02_voomLmFit_block.R")
echo "edgeRql:     $jid_eql ($PART_SMALL)"
echo "voomLmFit:   $jid_vlm ($PART_SMALL)"

# Arms 3 + 4 need per_study SE column AND are memory-heavy — bigmem.
jid_mfr=$(submit  metafor "$deps_step1" "$PART_HEAVY" 16 "$HEAVY_MEM" 48:00:00 "$SCRIPTS/03_metafor_REML_HKSJ.R")
jid_mash=$(submit mashr   "$deps_step1" "$PART_HEAVY" 16 "$HEAVY_MEM" 48:00:00 "$SCRIPTS/04_mashr_sharing.R" motifbreakr)
echo "metafor:     $jid_mfr ($PART_HEAVY $HEAVY_MEM)"
echo "mashr:       $jid_mash ($PART_HEAVY $HEAVY_MEM, motifbreakr env)"

# --- Tier 2: concordance after all 4 arms (mashr is afterany so failure doesn't block) ---
deps_conc="afterok:${jid_eql}:${jid_vlm}:${jid_mfr},afterany:${jid_mash}"
jid_conc=$(submit concord "$deps_conc" "$PART_SMALL" 4 32G 8:00:00 "$SCRIPTS/05_concordance.R")
echo "concordance: $jid_conc ($PART_SMALL)"

# --- Tier 3: figure + atlas refresh in parallel ---
jid_fig=$(submit  figS  "afterok:$jid_conc" "$PART_SMALL" 4 32G 4:00:00 "$PROJECT_ROOT/scripts/figures/figS_mega_validation.R")
jid_atlas=$(submit atlas "afterok:$jid_conc" "$PART_SMALL" 8 64G 4:00:00 "$PROJECT_ROOT/RNA-seq/27a_assemble_evidence_atlas.R")
echo "figS:        $jid_fig ($PART_SMALL)"
echo "atlas:       $jid_atlas ($PART_SMALL)"

echo
echo "All jobs submitted. Monitor with:"
echo "  squeue -u \$USER -n mvL_perstudy,mvL_eql,mvL_vlm,mvL_metafor,mvL_mashr,mvL_concord,mvL_figS,mvL_atlas"
