#!/bin/bash
#SBATCH --job-name=portal
#SBATCH --partition=cpu
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/portal/logs/rebuild_web_data_%j.log
# Big-RAM fallback (only if a step needs >210G): comment the `--partition=cpu`
# line above and uncomment the two lines below. bigmem+interactive enforces a
# ~500G QOSMinMemory floor, so also raise --mem to >=500G if you switch.
##SBATCH --partition=bigmem
##SBATCH --qos=interactive
# =============================================================================
# LEGACY PRE-RESOURCE web-data orchestrator. It emits rank-first convergence,
# progression, and Cas13-adjacent objects that are not part of the MASLD Gene Catalog
# contract. It fails closed unless explicitly invoked for provenance.
# =============================================================================
# Regenerates every artifact the public portal (masld-atlas-v2/) serves, from
# the current canonical results, in one re-runnable 7-phase DAG. Supersedes the
# partial scripts/run_atlas_rebuild_web.sbatch (27a/27b/preprocess only) and the
# retired scripts/figures/run_phase10_portal_refresh.sh.
#
# ENV note: `rnaseq` currently has a NumPy2/pandas ABI break, so R steps run
# under `rnaseq` and ALL Python steps under `spatial`. Each step names its own
# env via `micromamba run -n <env>` — there is no global `activate`.
#
# PARTITION note: defaults to cpu @96G (fits the cpu 210G cap; the natural home
# for this job and avoids the bigmem+interactive ~500G QOSMinMemory floor). A
# commented bigmem+interactive fallback is in the header for the rare big-RAM step.
#
# ENV-VAR CONTRACT (all optional; sensible defaults):
#   MASLD_PROJECT_ROOT  project root                  (default: this repo)
#   WEB_DATA_DIR        portal data output dir         (default: masld-atlas-v2/public/data)
#   HF_STAGE_DIR        HF Dataset staging dir         (default: <root>/masld-atlas-v2/hf_stage)
#   RUN_NETWORK         run network exporters (294-296) 1/0   (default: 1)
#   RUN_NEW_LAYERS      run new-layer generators        1/0   (default: 1)
#   PUSH                push to HF (never from job)     1/0   (default: 0)
#
# Idempotent: every generator overwrites its outputs; safe to re-run.
# Usage:  sbatch scripts/portal/rebuild_web_data.sbatch
# =============================================================================

set -eo pipefail   # NOT -u: micromamba/conda activation references unbound vars

if [[ "${ALLOW_LEGACY_PORTAL_REBUILD:-false}" != "true" ]]; then
  echo "REFUSED: legacy portal rebuild is outside the standalone Resource contract."
  echo "Use Plan 50/60; set ALLOW_LEGACY_PORTAL_REBUILD=true only for provenance-only regeneration."
  exit 64
fi

# ---------------------------------------------------------------------------
# Config / env-var contract
# ---------------------------------------------------------------------------
MASLD_PROJECT_ROOT="${MASLD_PROJECT_ROOT:-/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design}"
export MASLD_PROJECT_ROOT
BASE="${MASLD_PROJECT_ROOT}"
cd "${BASE}"

WEB_DATA_DIR="${WEB_DATA_DIR:-${BASE}/masld-atlas-v2/public/data}"
HF_STAGE_DIR="${HF_STAGE_DIR:-${BASE}/masld-atlas-v2/hf_dataset_stage}"
RUN_NETWORK="${RUN_NETWORK:-1}"
RUN_NEW_LAYERS="${RUN_NEW_LAYERS:-1}"
PUSH="${PUSH:-0}"

mkdir -p "${WEB_DATA_DIR}" "${BASE}/scripts/portal/logs"

# Convenience env-runners (steps name their own env; no global activate).
R_ENV="rnaseq"        # R only (27a/27b)
PY_ENV="spatial"      # all Python generators (rnaseq has a NumPy2/pandas ABI break)
run_R() {  micromamba run -n "${R_ENV}"  "$@"; }
run_py() { micromamba run -n "${PY_ENV}" "$@"; }

echo "=============================================================="
echo " MASLD web-data rebuild  —  $(date)  on $(hostname)"
echo "=============================================================="
echo "  MASLD_PROJECT_ROOT = ${MASLD_PROJECT_ROOT}"
echo "  WEB_DATA_DIR       = ${WEB_DATA_DIR}"
echo "  HF_STAGE_DIR       = ${HF_STAGE_DIR}"
echo "  RUN_NETWORK        = ${RUN_NETWORK}"
echo "  RUN_NEW_LAYERS     = ${RUN_NEW_LAYERS}"
echo "  PUSH               = ${PUSH}"
echo ""

# ---------------------------------------------------------------------------
# Step runner / result collector
# ---------------------------------------------------------------------------
succeeded=()
failed=()

run_step() {
  local name="$1"; shift
  echo "----------------------------------------------------------------"
  echo "[$(date '+%H:%M:%S')] STEP ${name}"
  echo "  cmd: $*"
  echo "----------------------------------------------------------------"
  if "$@"; then
    echo "[$(date '+%H:%M:%S')] ${name} OK"
    succeeded+=("${name}")
  else
    echo "[$(date '+%H:%M:%S')] ${name} FAILED (rc=$?)"
    failed+=("${name}")
  fi
  echo ""
}

# ---------------------------------------------------------------------------
# PREFLIGHT — fail fast on missing inputs / stale convergence
# ---------------------------------------------------------------------------
echo "=== PREFLIGHT ==="
ATLAS_CSV="${BASE}/RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
CONV_CSV="${BASE}/RNA-seq/results/multi_evidence/convergence_evidence.csv"
REQUIRED_INPUTS=(
  "${ATLAS_CSV}"
  "${CONV_CSV}"
  "${BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/canonical_deg_results.csv"
  "${BASE}/GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"
  "${BASE}/RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv"
)
preflight_ok=1
for f in "${REQUIRED_INPUTS[@]}"; do
  if [[ -e "$f" ]]; then
    echo "  OK   $f"
  else
    echo "  MISS $f"
    preflight_ok=0
  fi
done

# convergence_evidence.csv must be at least as fresh as the atlas it derives from
if [[ -e "${ATLAS_CSV}" && -e "${CONV_CSV}" ]]; then
  atlas_mtime=$(stat -c %Y "${ATLAS_CSV}")
  conv_mtime=$(stat -c %Y "${CONV_CSV}")
  if (( conv_mtime >= atlas_mtime )); then
    echo "  OK   convergence_evidence.csv mtime >= multi_evidence_atlas.csv mtime"
  else
    echo "  STALE convergence_evidence.csv is OLDER than multi_evidence_atlas.csv"
    echo "        (re-run RNA-seq convergence step before rebuilding web data)"
    preflight_ok=0
  fi
fi

if (( preflight_ok == 0 )); then
  echo "PREFLIGHT FAILED — aborting before any generator runs."
  exit 1
fi
echo "PREFLIGHT PASSED"
echo ""

# ===========================================================================
# PHASE 0 — Atlas rebuild + benchmark (R, env=rnaseq)
# ===========================================================================
echo "########## PHASE 0: atlas assembly (R) ##########"
run_step "27a_assemble_evidence_atlas" run_R Rscript "${BASE}/RNA-seq/27a_assemble_evidence_atlas.R"
run_step "27b_benchmark_presets"       run_R Rscript "${BASE}/RNA-seq/27b_benchmark_presets.R"

# ===========================================================================
# PHASE 1 — Stage-trajectory prerequisite CSV (Python, env=spatial)
# ===========================================================================
echo "########## PHASE 1: stage trajectories ##########"
run_step "build_stage_trajectories" run_py python "${BASE}/scripts/portal/build_stage_trajectories.py"

# ===========================================================================
# PHASE 2 — Core parquet + gene index + summary + featured (Python, spatial)
#   Produces: atlas.parquet, gene_index.json, atlas_summary.json, featured_genes.json
#   TODO(Team A): once the data contract lands, preprocess_atlas_data.py also
#                 emits the compact atlas_core.parquet (~40 UI columns).
# ===========================================================================
echo "########## PHASE 2: preprocess atlas -> parquet + json ##########"
run_step "preprocess_atlas_data" run_py python \
  "${BASE}/masld-atlas-v2/scripts/preprocess_atlas_data.py" --output-dir "${WEB_DATA_DIR}"

# ===========================================================================
# PHASE 3 — Phase3/4 + proteomics JSONs (Python, spatial; independent, ordered)
#   NOTE: generate_proteomics_data.py hardcodes its output to public/data
#         (does not honor WEB_DATA_DIR yet — parameterize later).
# ===========================================================================
echo "########## PHASE 3: phase3 / phase4 / proteomics ##########"
run_step "generate_phase3_data" run_py python \
  "${BASE}/masld-atlas-v2/scripts/generate_phase3_data.py" --output-dir "${WEB_DATA_DIR}"
run_step "generate_phase4_data" run_py python \
  "${BASE}/masld-atlas-v2/scripts/generate_phase4_data.py" --output-dir "${WEB_DATA_DIR}"
run_step "generate_proteomics_data" run_py python \
  "${BASE}/scripts/portal/generate_proteomics_data.py"
# drug_targets.parquet = full 27,943-gene drug universe (canonical pipeline
# headline 2/208/1504); separate from atlas_core's atlas-restricted column.
run_step "generate_drug_targets_data" run_py python \
  "${BASE}/scripts/portal/generate_drug_targets_data.py" --output-dir "${WEB_DATA_DIR}"

# ===========================================================================
# PHASE 4 — Convergence JSONs (Python, spatial)
#   MUST follow Phase 2: generate_convergence_data.py rewrites gene_index.json
#   in place (renames evidence keys, injects s5/s8), so gene_index.json from
#   preprocess must exist first.
# ===========================================================================
echo "########## PHASE 4: convergence (ranking + matrix) ##########"
run_step "generate_convergence_evidence_json" run_py python \
  "${BASE}/scripts/portal/generate_convergence_evidence_json.py"
run_step "generate_convergence_data" run_py python \
  "${BASE}/scripts/portal/generate_convergence_data.py"

# ===========================================================================
# PHASE 5 — NEW exploration-layer generators (Python, spatial)  [RUN_NEW_LAYERS]
#   All MASLD_PROJECT_ROOT-aware, take --output-dir, emit compact JSON + parquet.
#   Independent of each other; build_gene_detail_parquets replaces the retired
#   generate_gene_profiles.py (per-gene long-format parquets, DuckDB WHERE symbol=?).
# ===========================================================================
if [[ "${RUN_NEW_LAYERS}" == "1" ]]; then
  echo "########## PHASE 5: new layers ##########"
  run_step "generate_singlecell_data"     run_py python "${BASE}/scripts/portal/generate_singlecell_data.py"     --output-dir "${WEB_DATA_DIR}"
  run_step "generate_spatial_data"        run_py python "${BASE}/scripts/portal/generate_spatial_data.py"        --output-dir "${WEB_DATA_DIR}"
  run_step "generate_ancestry_coloc_data" run_py python "${BASE}/scripts/portal/generate_ancestry_coloc_data.py" --output-dir "${WEB_DATA_DIR}"
  run_step "generate_programs_data"       run_py python "${BASE}/scripts/portal/generate_programs_data.py"       --output-dir "${WEB_DATA_DIR}"
  run_step "generate_sex_data"            run_py python "${BASE}/scripts/portal/generate_sex_data.py"            --output-dir "${WEB_DATA_DIR}"
  run_step "build_gene_detail_parquets"   run_py python "${BASE}/scripts/portal/build_gene_detail_parquets.py"   --output-dir "${WEB_DATA_DIR}"
else
  echo "########## PHASE 5: new layers SKIPPED (RUN_NEW_LAYERS=0) ##########"
  echo ""
fi

# ===========================================================================
# PHASE 6 — Network export (Python, spatial)  [RUN_NETWORK]
#   Existing chain 294 -> 295 -> 296, then network_to_parquet.py which replaces
#   the 16,601 per-gene network JSONs with network_edges/nodes.parquet.
# ===========================================================================
if [[ "${RUN_NETWORK}" == "1" ]]; then
  echo "########## PHASE 6: network export (294 -> 295 -> 296 -> parquet) ##########"
  run_step "294_export_portal_v2"     run_py python "${BASE}/RNA-seq/294_export_portal_v2.py"
  run_step "295_build_filter_indexes" run_py python "${BASE}/RNA-seq/295_build_filter_indexes.py"
  run_step "296_community_subtype_labels" run_py python "${BASE}/RNA-seq/296_community_subtype_labels.py"
  run_step "network_to_parquet"       run_py python "${BASE}/scripts/portal/network_to_parquet.py" --output-dir "${WEB_DATA_DIR}"
else
  echo "########## PHASE 6: network SKIPPED (RUN_NETWORK=0) ##########"
  echo ""
fi

# ===========================================================================
# PHASE 7 — HF Dataset packaging (Python, stdlib)
#   stage_hf_dataset.py assembles HF_STAGE_DIR (loose parquet + compact JSON),
#   writes the dataset README card + .gitattributes, and PRINTS the `hf upload`
#   command. It NEVER pushes from the job (no compute-node egress); the user runs
#   the printed command from a host with egress. PUSH is advisory only here.
# ===========================================================================
echo "########## PHASE 7: HF packaging ##########"
run_step "stage_hf_dataset" run_py python "${BASE}/scripts/portal/stage_hf_dataset.py" \
  --web-dir "${WEB_DATA_DIR}" --stage-dir "${HF_STAGE_DIR}"
echo "  (PUSH=${PUSH}; the job only stages + prints the upload command — never pushes)"
echo ""

# ===========================================================================
# SUMMARY
# ===========================================================================
echo "=============================================================="
echo "=== SUMMARY  ($(date)) ==="
echo "=============================================================="
echo "Succeeded (${#succeeded[@]}):"
for s in "${succeeded[@]}"; do echo "  [ok]   $s"; done
if (( ${#failed[@]} > 0 )); then
  echo "FAILED (${#failed[@]}):"
  for f in "${failed[@]}"; do echo "  [FAIL] $f"; done
  echo "=== Completed WITH FAILURES ==="
  exit 1
fi
echo "=== All steps completed successfully ==="
