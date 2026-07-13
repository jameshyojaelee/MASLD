#!/bin/bash
# =============================================================================
# run_network_pipeline.sh — Bayesian Multiplex Gene Network Pipeline
# 4 Stages, 18 jobs, SLURM dependency chain
#
#   Stage 1 (R, cpu): 250 -> {251,252,253,254,255,256} in parallel
#   Stage 2 (Python, gpu+cpu): 260 -> 261 -> 262 -> 263 -> 264
#   Stage 3 (Python+R, cpu): 265 -> 266 -> {267, 268} in parallel
#   Stage 4 (Python+R, cpu): {267+268} -> 269; 268 -> 270
# =============================================================================

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
LOGDIR="${BASE}/logs/network_pipeline"
SCRIPTDIR="${LOGDIR}/scripts"
mkdir -p "${LOGDIR}" "${SCRIPTDIR}"

echo "=== Bayesian Multiplex Network Pipeline ==="
echo "Start: $(date)"
echo "Logs:  ${LOGDIR}/"
echo ""

submit_r() {
  local NUM=$1 NAME=$2 PART=$3 MEM=$4 CPUS=$5 SCRIPT=$6 DEP=${7:-}
  cat > "${SCRIPTDIR}/${NUM}.sbatch" << EOF
#!/bin/bash
#SBATCH --job-name=net_${NUM}_${NAME}
#SBATCH --partition=${PART}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=48:00:00

eval "\$(micromamba shell hook -s bash)"
micromamba activate rnaseq
cd ${BASE}
echo "${NUM} ${SCRIPT} | Start: \$(date)"
Rscript ${SCRIPT}
echo "${NUM} ${SCRIPT} | Exit: \$? | End: \$(date)"
EOF
  local DEP_FLAG=""
  [[ -n "$DEP" ]] && DEP_FLAG="--dependency=afterok:${DEP}"
  sbatch --parsable ${DEP_FLAG} \
    --output="${LOGDIR}/${NUM}_${NAME}_%j.out" \
    --error="${LOGDIR}/${NUM}_${NAME}_%j.err" \
    "${SCRIPTDIR}/${NUM}.sbatch"
}

submit_py() {
  local NUM=$1 NAME=$2 PART=$3 MEM=$4 CPUS=$5 SCRIPT=$6 DEP=${7:-} GPU=${8:-}
  local GPU_LINE=""
  [[ -n "$GPU" ]] && GPU_LINE="#SBATCH --gres=gpu:1"
  cat > "${SCRIPTDIR}/${NUM}.sbatch" << EOF
#!/bin/bash
#SBATCH --job-name=net_${NUM}_${NAME}
#SBATCH --partition=${PART}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=48:00:00
${GPU_LINE}

eval "\$(micromamba shell hook -s bash)"
micromamba activate spatial
cd ${BASE}
echo "${NUM} ${SCRIPT} | Start: \$(date)"
python ${SCRIPT}
echo "${NUM} ${SCRIPT} | Exit: \$? | End: \$(date)"
EOF
  local DEP_FLAG=""
  [[ -n "$DEP" ]] && DEP_FLAG="--dependency=afterok:${DEP}"
  sbatch --parsable ${DEP_FLAG} \
    --output="${LOGDIR}/${NUM}_${NAME}_%j.out" \
    --error="${LOGDIR}/${NUM}_${NAME}_%j.err" \
    "${SCRIPTDIR}/${NUM}.sbatch"
}

# ── Stage 1: Node definition + Edge construction ──
J250=$(submit_r 250 nodes cpu 16G 4 250_define_network_nodes.R)
echo "250 nodes:       ${J250}"

J251=$(submit_r 251 ppi      cpu 32G 4  251_ppi_edges.R                   "${J250}")
J252=$(submit_r 252 coexpr   cpu 64G 16 252_coexpression_edges.R          "${J250}")
J253=$(submit_r 253 reg_lr   cpu 16G 4  253_regulon_lr_cerna_edges.R      "${J250}")
J254=$(submit_r 254 genetic  cpu 32G 4  254_genetic_coloc_edges.R         "${J250}")
J255=$(submit_r 255 pathway  cpu 64G 8  255_pathway_edges.R               "${J250}")
J256=$(submit_r 256 spat_cos cpu 32G 8  256_spatial_cosmos_xspecies_edges.R "${J250}")
echo "251 ppi:         ${J251}  (after 250)"
echo "252 coexpr:      ${J252}  (after 250)"
echo "253 reg/lr/cerna: ${J253}  (after 250)"
echo "254 genetic:     ${J254}  (after 250)"
echo "255 pathway:     ${J255}  (after 250)"
echo "256 spat/cos/xsp: ${J256}  (after 250)"

# ── Stage 2: Bayesian calibration ──
EDGE_DEPS="${J251}:${J252}:${J253}:${J254}:${J255}:${J256}"
J260=$(submit_py 260 permnull gpu 64G 8 260_permutation_nulls.py       "${EDGE_DEPS}" gpu)
# 261 (KDE null density) + 262 (KDE two-group posteriors) RETIRED 2026-07-04
# (round-2 audit B4d): the KDE-lfdr posterior was degenerate (pi0->0.99,
# posteriors->0; see posterior_summary.csv) and is superseded by 262c direct
# confidence scores. Scripts + degenerate outputs archived under
# data/archive/kde_lfdr_retired_2026-07-04/. 262c reads the raw stage-1 edges
# and writes the same posterior_edges/ filenames 263 consumes.
J262=$(submit_py 262c confid   cpu 32G 4 262c_confidence_scores.py      "${J260}")
J263=$(submit_py 263 noisyor  cpu 32G 8 263_noisy_or_composite.py       "${J262}")
J264=$(submit_py 264 diagnost cpu 16G 4 264_bayesian_diagnostics.py     "${J263}")
echo "260 permnull:    ${J260}  (after 251-256)"
echo "262c confidence: ${J262}  (after 260; replaces retired 261+262 KDE-lfdr)"
echo "263 noisy-or:    ${J263}  (after 262c)"
echo "264 diagnostics: ${J264}  (after 263)"

# ── Stage 3: Graph assembly ──
J265=$(submit_py 265 assemble cpu 64G 8 265_assemble_multiplex_graph.py "${J263}")
J266=$(submit_py 266 leiden   cpu 32G 8 266_hierarchical_leiden.py      "${J265}")
J267=$(submit_py 267 neighbor cpu 32G 8 267_gene_neighborhoods.py       "${J266}")
J268=$(submit_r  268 enrich   cpu 32G 8 268_community_enrichment.R      "${J266}")
echo "265 assemble:    ${J265}  (after 263)"
echo "266 leiden:      ${J266}  (after 265)"
echo "267 neighborhoods: ${J267}  (after 266)"
echo "268 enrichment:  ${J268}  (after 266)"

# ── Stage 4: Export ──
J269=$(submit_py 269 export   cpu 32G 8 269_export_portal_json.py       "${J267}:${J268}")
J270=$(submit_r  270 atlas    cpu 16G 4 270_atlas_integration.R         "${J268}")
echo "269 export:      ${J269}  (after 267+268)"
echo "270 atlas:       ${J270}  (after 268)"

echo ""
echo "=== Pipeline Submitted (18 jobs) ==="
echo "Monitor: squeue -u \$USER | grep net_"
echo "Logs:    ${LOGDIR}/"
echo "End:     $(date)"
