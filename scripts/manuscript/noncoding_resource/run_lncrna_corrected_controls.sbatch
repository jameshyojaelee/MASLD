#!/bin/bash
#SBATCH --job-name=nc_lnc_corr
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nc_lnc_corr_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/nc_lnc_corr_%j.err

# Figure 3 lncRNA inset on the corrected-control F_five fit (review item 1, Phase 3).
# The GENCODE v49 identity, lncRNA classes and strandedness source gate do not
# depend on control labels and are reused. The two producers rerun on the
# corrected arm, then fig3_lncrna_role.R renders from them:
#   1 build_lncrna_biotype_association.R  (about 10 s)
#   2 build_bulk_lncrna_candidate.R       (per-cohort, leave-one-cohort-out and
#                                          equal-weight refits; 1 h 35 min, 2.4 GB
#                                          in job 19810556)
#   3 fig3_lncrna_role.R                  (seconds)
#
#   sbatch --dependency=afterok:<fcorr_repro job> --export=ALL,\
# FIT_ARM=<RNA-seq>/results/manuscript_release/candidates/bulk-corrected-controls-<ts>/primary_source_controls,\
# PRIMARY_ARM=<RNA-seq>/results/manuscript_release/candidates/bulk-corrected-controls-<ts>-reproduction/primary_source_controls,\
# LNCRNA_ROOT=<RNA-seq>/results/noncoding_resource/candidates/noncoding-bulk-corrected-controls-<UTC ts>,\
# FIG3_LNCRNA_OUT=<BASE>/figures/candidates/fig3-lncrna-corrected-controls-<UTC ts> \
#     scripts/manuscript/noncoding_resource/run_lncrna_corrected_controls.sbatch
set -euo pipefail
umask 027
: "${FIT_ARM:?Set the corrected fit-root arm (merged_dge.rds, meta_matched.rds)}"
: "${PRIMARY_ARM:?Set the corrected reproduction-root arm (deg_results.csv, VALIDATED.json)}"
: "${LNCRNA_ROOT:?Set a new noncoding candidate root}"
: "${FIG3_LNCRNA_OUT:?Set a new figure candidate root}"
test ! -e "${LNCRNA_ROOT}"
test ! -e "${FIG3_LNCRNA_OUT}"
for f in "${FIT_ARM}/merged_dge.rds" "${FIT_ARM}/meta_matched.rds" "${PRIMARY_ARM}/deg_results.csv" \
         "${PRIMARY_ARM}/VALIDATED.json" "${PRIMARY_ARM}/model_design.tsv" "${PRIMARY_ARM}/expected_counts.json"; do
  test -s "$f" || { echo "Missing input: $f" >&2; exit 2; }
done

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
RNASEQ_ENV="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
NC="${PROJECT_ROOT}/RNA-seq/results/noncoding_resource/candidates/noncoding-bulk-2026-08-11"
IDENTITY_ROOT="${NC}/annotation/GENCODE_v49/run_19749483"
SOURCE_GATE="${NC}/source_gate/run_19751171/source_gate_status.tsv"
ADOPTED_BIOTYPE="${NC}/biotype_association/run_canonical_gate_2026-08-12"
ADOPTED_BULK="${NC}/bulk_disease/run_canonical_gate_2026-08-12"
SRC="${PROJECT_ROOT}/scripts/manuscript/noncoding_resource"
export MASLD_PROJECT_ROOT="${PROJECT_ROOT}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export PATH="${RNASEQ_ENV}/bin:/gpfs/commons/home/jameslee/.local/bin:${PATH}"
export CONDA_PREFIX="${RNASEQ_ENV}"
export PYTHONNOUSERSITE=1
mkdir "${LNCRNA_ROOT}"

"${RNASEQ_ENV}/bin/Rscript" --vanilla "${SRC}/build_lncrna_biotype_association.R" \
  "--source-gate=${SOURCE_GATE}" \
  "--gene-identity=${IDENTITY_ROOT}/gene_identity.tsv" \
  "--fit-arm=${FIT_ARM}" "--primary-arm=${PRIMARY_ARM}" \
  "--output=${LNCRNA_ROOT}/biotype_association"

"${RNASEQ_ENV}/bin/Rscript" --vanilla "${SRC}/build_bulk_lncrna_candidate.R" \
  "--source-gate=${SOURCE_GATE}" \
  "--gene-identity=${IDENTITY_ROOT}/gene_identity.tsv" \
  "--lncrna-class=${IDENTITY_ROOT}/lncrna_genomic_class.tsv" \
  "--fit-arm=${FIT_ARM}" "--primary-arm=${PRIMARY_ARM}" \
  "--output=${LNCRNA_ROOT}/bulk_disease"

FIG3_LNCRNA_BIOTYPE_DIR="${LNCRNA_ROOT}/biotype_association" \
FIG3_LNCRNA_BULK_DIR="${LNCRNA_ROOT}/bulk_disease" \
FIG3_LNCRNA_EXPECT_ADOPTED=0 \
FIG3_LNCRNA_OUT="${FIG3_LNCRNA_OUT}" \
  "${RNASEQ_ENV}/bin/Rscript" --vanilla "${PROJECT_ROOT}/scripts/figures/fig3_lncrna_role.R"

PDF="${FIG3_LNCRNA_OUT}/panels/fig3_lncrna_role_candidate.pdf"
test "$(pdfinfo "${PDF}" | awk '/^Pages:/ {print $2}')" = "1"
! strings "${PDF}" | grep -q "/Subtype /Type3"
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${PDF}"
mkdir -p "${FIG3_LNCRNA_OUT}/proofs"
pdftotext "${PDF}" "${FIG3_LNCRNA_OUT}/proofs/fig3_lncrna_role_candidate.txt"
pdftoppm -png -r 400 -singlefile "${PDF}" "${FIG3_LNCRNA_OUT}/proofs/fig3_lncrna_role_candidate"

# Adopted versus corrected: biotype census, the association, and every summary value.
"${RNASEQ_ENV}/bin/python3" - "${ADOPTED_BIOTYPE}" "${LNCRNA_ROOT}/biotype_association" \
  "${ADOPTED_BULK}" "${LNCRNA_ROOT}/bulk_disease" \
  "${PROJECT_ROOT}/figures/main/fig3_bulk_transcriptomics/source_tables/fig3_lncrna_role/summary.tsv" \
  "${FIG3_LNCRNA_OUT}/source_tables/fig3_lncrna_role/summary.tsv" \
  "${FIG3_LNCRNA_OUT}/old_vs_new.tsv" <<'PY'
import sys
import pandas as pd
old_bt, new_bt, old_bulk, new_bulk, old_sum, new_sum, out = sys.argv[1:]
rows = []
for name in ("bulk_deg_counts_by_biotype.tsv",):
    o, n = (pd.read_csv(f"{d}/{name}", sep="\t").set_index("display_biotype") for d in (old_bt, new_bt))
    for b in o.index:
        for c in o.columns:
            rows.append({"table": name, "row": b, "field": c, "adopted": o.at[b, c], "corrected": n.at[b, c]})
o, n = (pd.read_csv(f"{d}/biotype_deg_association.tsv", sep="\t").iloc[0] for d in (old_bt, new_bt))
for c in ("n_genes", "n_lncrna", "n_protein_coding", "odds_ratio", "ci_lower", "ci_upper", "nominal_p"):
    rows.append({"table": "biotype_deg_association.tsv", "row": "lncRNA_vs_protein_coding",
                 "field": c, "adopted": o[c], "corrected": n[c]})
o, n = (pd.read_csv(f"{d}/bulk_lncrna_verdict.tsv", sep="\t").set_index("metric")["value"]
        for d in (old_bulk, new_bulk))
for m in o.index:
    rows.append({"table": "bulk_lncrna_verdict.tsv", "row": m, "field": "value",
                 "adopted": o[m], "corrected": n[m]})
o, n = (pd.read_csv(p, sep="\t") for p in (old_sum, new_sum))
key = ["section", "metric"]
m = o.merge(n, on=key, how="outer", suffixes=("_adopted", "_corrected"))
for _, r in m.iterrows():
    for c in ("value", "lower", "upper", "n"):
        rows.append({"table": "summary.tsv", "row": f"{r['section']}:{r['metric']}", "field": c,
                     "adopted": r[f"{c}_adopted"], "corrected": r[f"{c}_corrected"]})
table = pd.DataFrame(rows)
table.to_csv(out, sep="\t", index=False)
print(table.to_string(index=False))
PY

python3 "${PROJECT_ROOT}/scripts/manuscript/validate_resource_scope.py"
echo "LNCRNA_CORRECTED_COMPLETE ${LNCRNA_ROOT} ${FIG3_LNCRNA_OUT}"
