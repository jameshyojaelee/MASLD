#!/bin/bash
set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq"
LOGDIR="${BASE}/logs/network_pipeline"
SCRIPTDIR="${LOGDIR}/scripts"
mkdir -p "${LOGDIR}" "${SCRIPTDIR}"

echo "=== Stages 2-4 (Stage 1 already complete) ==="
echo "Start: $(date)"

write_sbatch() {
  local NUM=$1 NAME=$2 PART=$3 MEM=$4 CPUS=$5 ENV=$6 CMD=$7
  cat > "${SCRIPTDIR}/${NUM}_${NAME}.sbatch" << EOF
#!/bin/bash
#SBATCH --job-name=net_${NUM}_${NAME}
#SBATCH --partition=${PART}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=48:00:00

eval "\$(micromamba shell hook -s bash)"
micromamba activate ${ENV}
cd ${BASE}
echo "${NUM} | Start: \$(date)"
${CMD}
echo "${NUM} | Exit: \$? | End: \$(date)"
EOF
}

write_sbatch 260 permnull  cpu 64G 8 spatial "python 260_permutation_nulls.py"
write_sbatch 261 nulldens  cpu 16G 4 spatial "python 261_null_density_estimation.py"
write_sbatch 262 postcomp  cpu 32G 4 spatial "python 262_compute_posteriors.py"
write_sbatch 263 noisyor   cpu 32G 8 spatial "python 263_noisy_or_composite.py"
write_sbatch 264 diagnost  cpu 16G 4 spatial "python 264_bayesian_diagnostics.py"
write_sbatch 265 assemble  cpu 64G 8 spatial "python 265_assemble_multiplex_graph.py"
write_sbatch 266 leiden    cpu 32G 8 spatial "python 266_hierarchical_leiden.py"
write_sbatch 267 neighbor  cpu 32G 8 spatial "python 267_gene_neighborhoods.py"
write_sbatch 268 enrich    cpu 32G 8 rnaseq  "Rscript 268_community_enrichment.R"
write_sbatch 269 export    cpu 32G 8 spatial "python 269_export_portal_json.py"
write_sbatch 270 atlas     cpu 16G 4 rnaseq  "Rscript 270_atlas_integration.R"

J260=$(sbatch --parsable --output="${LOGDIR}/260_%j.out" --error="${LOGDIR}/260_%j.err" "${SCRIPTDIR}/260_permnull.sbatch")
echo "260 permnull:    ${J260}"

J261=$(sbatch --parsable --dependency=afterok:${J260} --output="${LOGDIR}/261_%j.out" --error="${LOGDIR}/261_%j.err" "${SCRIPTDIR}/261_nulldens.sbatch")
J262=$(sbatch --parsable --dependency=afterok:${J261} --output="${LOGDIR}/262_%j.out" --error="${LOGDIR}/262_%j.err" "${SCRIPTDIR}/262_postcomp.sbatch")
J263=$(sbatch --parsable --dependency=afterok:${J262} --output="${LOGDIR}/263_%j.out" --error="${LOGDIR}/263_%j.err" "${SCRIPTDIR}/263_noisyor.sbatch")
J264=$(sbatch --parsable --dependency=afterok:${J263} --output="${LOGDIR}/264_%j.out" --error="${LOGDIR}/264_%j.err" "${SCRIPTDIR}/264_diagnost.sbatch")
J265=$(sbatch --parsable --dependency=afterok:${J263} --output="${LOGDIR}/265_%j.out" --error="${LOGDIR}/265_%j.err" "${SCRIPTDIR}/265_assemble.sbatch")
J266=$(sbatch --parsable --dependency=afterok:${J265} --output="${LOGDIR}/266_%j.out" --error="${LOGDIR}/266_%j.err" "${SCRIPTDIR}/266_leiden.sbatch")
J267=$(sbatch --parsable --dependency=afterok:${J266} --output="${LOGDIR}/267_%j.out" --error="${LOGDIR}/267_%j.err" "${SCRIPTDIR}/267_neighbor.sbatch")
J268=$(sbatch --parsable --dependency=afterok:${J266} --output="${LOGDIR}/268_%j.out" --error="${LOGDIR}/268_%j.err" "${SCRIPTDIR}/268_enrich.sbatch")
J269=$(sbatch --parsable --dependency=afterok:${J267}:${J268} --output="${LOGDIR}/269_%j.out" --error="${LOGDIR}/269_%j.err" "${SCRIPTDIR}/269_export.sbatch")
J270=$(sbatch --parsable --dependency=afterok:${J268} --output="${LOGDIR}/270_%j.out" --error="${LOGDIR}/270_%j.err" "${SCRIPTDIR}/270_atlas.sbatch")

echo "261 nulldens:    ${J261}  (after 260)"
echo "262 posteriors:  ${J262}  (after 261)"
echo "263 noisy-or:    ${J263}  (after 262)"
echo "264 diagnostics: ${J264}  (after 263)"
echo "265 assemble:    ${J265}  (after 263)"
echo "266 leiden:      ${J266}  (after 265)"
echo "267 neighborhoods: ${J267}  (after 266)"
echo "268 enrichment:  ${J268}  (after 266)"
echo "269 export:      ${J269}  (after 267+268)"
echo "270 atlas:       ${J270}  (after 268)"
echo ""
echo "=== Submitted 11 jobs ==="
echo "End: $(date)"
