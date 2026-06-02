#!/bin/bash -l
# refetch_okur2026_supplement.sh
# Automated re-fetch of the Okur et al. 2026 supplementary data from Zenodo
# once the embargo lifts on 2026-05-01 (DOI: 10.5281/zenodo.18662736).
#
# Steps:
#   1. Check Zenodo API for embargo status.
#   2. If open, download the full supplement archive.
#   3. Extract Supplemental Table 1 (full 90-gene panel).
#   4. Replace data/external/okur2026_monogenic/okur2026_panel_genes_partial.tsv
#      with the full list as okur2026_panel_genes.tsv.
#   5. Re-run RNA-seq/53_monogenic_convergence.R and figS_monogenic_convergence.R.
#
# Schedule: cron or manual run after 2026-05-01.
# Example cron entry (daily at 09:00 until successful):
#   0 9 * * * /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/scripts/refetch_okur2026_supplement.sh >> /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/data/external/okur2026_monogenic/refetch.log 2>&1

set -o pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

BASE=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT_DIR="${BASE}/data/external/okur2026_monogenic"
DONE_FLAG="${OUT_DIR}/full_panel_fetched.flag"
cd "${BASE}"

if [ -f "${DONE_FLAG}" ]; then
  echo "[$(date)] Full panel already fetched and processed. Exiting."
  exit 0
fi

ZENODO_API="https://zenodo.org/api/records/18662736"

echo "[$(date)] Probing Zenodo record..."
META=$(curl -sL --max-time 30 "${ZENODO_API}")
EMBARGO=$(echo "${META}" | python3 -c "import sys,json; d=json.load(sys.stdin); m=d.get('metadata',{}); print(m.get('access_right','unknown'),'|',m.get('embargo_date','-'))")
echo "[$(date)] Zenodo status: ${EMBARGO}"

if echo "${META}" | grep -q '"access_right": *"embargoed"'; then
  echo "[$(date)] Still embargoed. Will retry later."
  exit 0
fi

echo "[$(date)] Embargo lifted. Downloading archive..."
ARCHIVE="${OUT_DIR}/supplement.zip"
curl -sL --max-time 600 -o "${ARCHIVE}" "https://zenodo.org/api/records/18662736/files-archive"
if [ ! -s "${ARCHIVE}" ]; then
  echo "[$(date)] ERROR: archive download empty" >&2
  exit 1
fi

echo "[$(date)] Extracting..."
mkdir -p "${OUT_DIR}/zenodo_extract"
(cd "${OUT_DIR}/zenodo_extract" && unzip -o "${ARCHIVE}")
ls "${OUT_DIR}/zenodo_extract/"

# Heuristic: look for "Supplemental Table 1" or "Panel" xlsx/csv/tsv
echo "[$(date)] Searching for panel file..."
PANEL_FILE=$(find "${OUT_DIR}/zenodo_extract" -type f \( -iname "*table*1*" -o -iname "*panel*" -o -iname "*gene*" \) \
             \( -iname "*.xlsx" -o -iname "*.csv" -o -iname "*.tsv" \) | head -1)
if [ -z "${PANEL_FILE}" ]; then
  echo "[$(date)] WARNING: could not auto-locate panel file in archive. Manual review needed."
  ls -la "${OUT_DIR}/zenodo_extract/"
  exit 2
fi
echo "[$(date)] Panel file: ${PANEL_FILE}"

# Convert xlsx → tsv if needed (using R + readxl)
if echo "${PANEL_FILE}" | grep -qi "\.xlsx$"; then
  PANEL_TSV="${OUT_DIR}/okur2026_panel_genes.tsv"
  Rscript -e "d <- readxl::read_xlsx('${PANEL_FILE}', sheet=1); write.table(d, '${PANEL_TSV}', sep='\t', quote=FALSE, row.names=FALSE)"
else
  cp "${PANEL_FILE}" "${OUT_DIR}/okur2026_panel_genes.tsv"
fi
echo "[$(date)] Full panel: $(wc -l < ${OUT_DIR}/okur2026_panel_genes.tsv) lines"

# Archive the partial file
if [ -f "${OUT_DIR}/okur2026_panel_genes_partial.tsv" ] && [ ! -f "${OUT_DIR}/okur2026_panel_genes_partial.tsv.pre_fullfetch" ]; then
  mv "${OUT_DIR}/okur2026_panel_genes_partial.tsv" \
     "${OUT_DIR}/okur2026_panel_genes_partial.tsv.pre_fullfetch"
fi

# Re-run convergence scripts
echo "[$(date)] Re-running 53_monogenic_convergence.R + figS_monogenic_convergence.R..."
Rscript RNA-seq/53_monogenic_convergence.R
Rscript scripts/figures/figS_monogenic_convergence.R

touch "${DONE_FLAG}"
echo "[$(date)] Done. Atlas refreshed with full 90-gene panel."
