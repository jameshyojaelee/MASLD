#!/bin/bash
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
# Download the 224 C1-matched GTEx liver slides (IDC DICOM) to data/external/gtex_liver_slides/series/<SeriesInstanceUID>/.
# Usage: sbatch --job-name=<n> --output=<log>/%j.out c2_download_slides.sbatch <c1_exec_dir>
set -euo pipefail
repo=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
dest=$repo/data/external/gtex_liver_slides/series
mkdir -p "$dest"
export PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/.local/bin/micromamba run -n he_molecular python - "$1/liver_slides_matched.tsv" "$dest" <<'PY'
import sys, pandas as pd
from idc_index import IDCClient
m = pd.read_csv(sys.argv[1], sep="\t")
c = IDCClient()
c.download_from_selection(seriesInstanceUID=list(m["SeriesInstanceUID"]), downloadDir=sys.argv[2], dirTemplate="%SeriesInstanceUID")
PY
find "$dest" -maxdepth 1 -mindepth 1 -type d | wc -l > "$dest/../series_count.txt"
du -s --apparent-size "$dest" >> "$dest/../series_count.txt"
echo slides done
