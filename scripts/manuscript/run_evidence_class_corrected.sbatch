#!/bin/bash
#SBATCH --job-name=evclass_corr
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/evidence_class/evclass_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/evidence_class/evclass_%j.err

# Evidence classes on the corrected bulk fit and the corrected COLOC release.
# Required:
#   EVIDENCE_CLASS_DEG     corrected F_five deg_results.csv
#   EVIDENCE_CLASS_COLOC   susie_coloc_r2_<ts>/susie_coloc_all_gwas.csv
#   EVIDENCE_CLASS_OUT_DIR new directory; the manifest and claim ledger go here,
#                          not to docs/manuscript/release
#   MANUSCRIPT_RELEASE_ID  label written into every output row
# The atlas and convergence inputs keep their defaults: the builder reads only
# gene identity, biotype, and the proteomic/spatial/single-cell endpoints from
# the atlas, and only hashes convergence_evidence.csv.
set -euo pipefail
umask 027

project_root="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
rnaseq_env="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq"
: "${EVIDENCE_CLASS_DEG:?EVIDENCE_CLASS_DEG is required}"
: "${EVIDENCE_CLASS_COLOC:?EVIDENCE_CLASS_COLOC is required}"
: "${EVIDENCE_CLASS_OUT_DIR:?EVIDENCE_CLASS_OUT_DIR is required}"
: "${MANUSCRIPT_RELEASE_ID:?MANUSCRIPT_RELEASE_ID is required}"
test -s "${EVIDENCE_CLASS_DEG}"
test -s "${EVIDENCE_CLASS_COLOC}"
test ! -e "${EVIDENCE_CLASS_OUT_DIR}"

export MASLD_PROJECT_ROOT="${project_root}"
export EVIDENCE_CLASS_DEG EVIDENCE_CLASS_COLOC EVIDENCE_CLASS_OUT_DIR MANUSCRIPT_RELEASE_ID
export LANG="en_US.UTF-8"
export LC_ALL="en_US.UTF-8"

"${rnaseq_env}/bin/Rscript" --vanilla \
  "${project_root}/scripts/manuscript/build_evidence_class_release.R"
sha256sum "${EVIDENCE_CLASS_DEG}" "${EVIDENCE_CLASS_COLOC}" > "${EVIDENCE_CLASS_OUT_DIR}/input_sha256.txt"
echo "EVCLASS_COMPLETE ${EVIDENCE_CLASS_OUT_DIR}"
