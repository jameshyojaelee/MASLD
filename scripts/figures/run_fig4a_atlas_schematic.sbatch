#!/usr/bin/env bash
#SBATCH --job-name=fig4a_schematic
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=logs/fig4a_schematic_%j.out
#SBATCH --error=logs/fig4a_schematic_%j.err
set -euo pipefail
umask 027
project_root=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
run_stamp=$(date -u +%Y%m%dT%H%M%SZ)
out="${project_root}/figures/candidates/fig4a-atlas-schematic-${run_stamp}_${SLURM_JOB_ID}"
rscript_bin=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq_r44/bin/Rscript
python_bin=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
pdf_tools=/gpfs/commons/home/jameslee/micromamba/bin
test ! -e "${out}"
mkdir -p "${out}/snapshot"
cp -p scripts/figures/fig4a_atlas_program_schematic.R scripts/figures/run_fig4a_atlas_schematic.sbatch \
    scripts/figures/normalize_pdf_page_box.py docs/RESULTS.md "${out}/snapshot/"
printf '%s\n' "${out}" > "logs/fig4a_schematic_${SLURM_JOB_ID}.root"
"${rscript_bin}" "${out}/snapshot/fig4a_atlas_program_schematic.R" "${project_root}" "${out}"
pdf="${out}/fig4a_atlas_program_schematic.pdf"
"${python_bin}" -B "${out}/snapshot/normalize_pdf_page_box.py" "${pdf}" 2.5 0.4027777777777777777777777778
"${pdf_tools}/pdfinfo" "${pdf}" > "${out}/pdfinfo.txt"
"${pdf_tools}/pdftotext" "${pdf}" "${out}/artwork_text.txt"
"${pdf_tools}/pdffonts" "${pdf}" > "${out}/fonts.txt"
gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
"${python_bin}" -B - "${out}" <<'PY'
from pathlib import Path
import re, sys
out = Path(sys.argv[1])
info = (out / 'pdfinfo.txt').read_text()
assert re.search(r'^Pages:\s+1\s*$', info, re.M)
size = re.search(r'^Page size:\s+([\d.]+) x ([\d.]+) pts', info, re.M)
assert size and tuple(map(float, size.groups())) == (180.0, 29.0)
text = (out / 'artwork_text.txt').read_text()
assert all(value in text for value in ('7 datasets', '104 donors', '1.23M cells', '102 analyzed donors', '117 programs', '5 lineages'))
fonts = (out / 'fonts.txt').read_text()
assert 'Helvetica' in fonts and 'Type 3' not in fonts
(out / 'COMPLETE.txt').write_text('Figure 4A schematic source, page, text, font and PDF checks passed.\n')
print('FIG4A_CHECK_PASS')
PY
sha256sum "${out}"/*.pdf "${out}"/*.tsv "${out}/snapshot/"* > "${out}/checksums.txt"
