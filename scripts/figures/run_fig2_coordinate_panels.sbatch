#!/usr/bin/env bash
# Render the adopted coordinate-corrected Figure 2 estimates as individual PDFs.
#SBATCH --job-name=fig2_coord_pdf
#SBATCH --account=nslab
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/fig2_coordinate_panels_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/fig2_coordinate_panels_%j.err

set -euo pipefail
umask 027
project_root=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
analysis_root="${project_root}/figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816"
python_bin=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
rscript_bin=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq_r44/bin/Rscript
pdf_tools_bin=/gpfs/commons/home/jameslee/micromamba/bin
export PATH="${pdf_tools_bin}:${PATH}"
test -x "${pdf_tools_bin}/pdffonts"
run_stamp=$(date -u +%Y%m%dT%H%M%SZ)
run_root="${project_root}/figures/candidates/fig2-coordinate-panels-${run_stamp}_${SLURM_JOB_ID}"
test -s "${analysis_root}/COMPLETE.txt"
test ! -e "${run_root}"
mkdir "${run_root}"
mkdir -p "${run_root}/snapshot/scripts/figures" "${run_root}/snapshot/trait_directness" "${run_root}/snapshot/multiplicity"
printf '%s\n' "${run_root}" > "${project_root}/logs/fig2_coordinate_panels_${SLURM_JOB_ID}.run_root.txt"
printf '[fig2] job=%s host=%s run=%s start=%s\n' "${SLURM_JOB_ID}" "$(hostname)" "${run_root}" "$(date -u)"

for script_name in fig2_trait_directness_panels.R fig2_consequence_palette.R \
    normalize_pdf_page_box.py validate_fig2d_pip_architecture.py; do
    cp -p "${project_root}/scripts/figures/${script_name}" "${run_root}/snapshot/scripts/figures/"
done
cp -p "${project_root}/scripts/figures/run_fig2_coordinate_panels.sbatch" "${run_root}/launcher.sbatch"
for source_name in pip_architecture_estimates.tsv trait_class_contrasts.tsv \
    portfolio_independence.tsv direct_trait_coding_concentration.tsv \
    input_manifest.tsv environment.txt; do
    cp -p "${analysis_root}/trait_directness/${source_name}" "${run_root}/snapshot/trait_directness/"
done
cp -p "${analysis_root}/multiplicity/fig2_family_multiplicity.tsv" "${run_root}/snapshot/multiplicity/"

"${python_bin}" -B - "${project_root}" "${analysis_root}" "${run_root}" <<'PY'
import csv
import hashlib
from pathlib import Path
import sys

root, analysis, run = map(Path, sys.argv[1:])
paths = [root / 'scripts/figures' / name for name in (
    'fig2_trait_directness_panels.R', 'fig2_consequence_palette.R',
    'normalize_pdf_page_box.py', 'validate_fig2d_pip_architecture.py',
    'run_fig2_coordinate_panels.sbatch')]
paths.extend(analysis / 'trait_directness' / name for name in (
    'pip_architecture_estimates.tsv', 'trait_class_contrasts.tsv',
    'portfolio_independence.tsv', 'direct_trait_coding_concentration.tsv',
    'input_manifest.tsv', 'environment.txt'))
paths.append(analysis / 'multiplicity/fig2_family_multiplicity.tsv')
paths.extend(path for path in (run / 'snapshot').rglob('*') if path.is_file())
with (run / 'input_checksums_before.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=('path', 'sha256'), delimiter='\t')
    writer.writeheader()
    for path in sorted(set(paths)):
        writer.writerow({'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
PY

export MASLD_PROJECT_ROOT="${project_root}"
export MASLD_FIGURE_PYTHON="${python_bin}"
export FIG2_TRAIT_DIRECTNESS_STATS="${run_root}/snapshot/trait_directness"
export FIG2_CORRECTED_INFERENCE_DIR="${run_root}/snapshot/multiplicity"
export FIG2_TRAIT_DIRECTNESS_OUT="${run_root}"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"${rscript_bin}" -e 'stopifnot(capabilities("cairo")); library(data.table); sessionInfo()' \
    > "${run_root}/render_environment.txt" 2>&1
"${rscript_bin}" "${run_root}/snapshot/scripts/figures/fig2_trait_directness_panels.R" \
    > "${run_root}/render.log" 2>&1

for pdf in "${run_root}"/*.pdf; do
    stem="$(basename "${pdf}" .pdf)"
    pdfinfo "${pdf}" > "${run_root}/${stem}_pdfinfo.txt"
    pdftotext "${pdf}" "${run_root}/${stem}_artwork_text.txt"
    pdffonts "${pdf}" > "${run_root}/${stem}_fonts.txt"
    gs -q -dNOPAUSE -dBATCH -sDEVICE=nullpage "${pdf}"
done

"${python_bin}" -B - "${run_root}" "${SLURM_JOB_ID}" <<'PY'
import csv
import hashlib
from pathlib import Path
import re
import sys

run = Path(sys.argv[1])
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
stem = 'Fig2D_pip_architecture_by_trait_directness'
fields = {
    'panel_index': 'Figure 2D', 'promotion_date': '2026-10-13',
    'requested_candidate': str(run / f'{stem}.pdf'),
    'promoted_artwork_source': str(run / f'{stem}.pdf'),
    'promoted_pdf_sha256': digest(run / f'{stem}.pdf'),
    'source_estimates_sha256': digest(run / f'{stem}_source.tsv'),
    'source_contrasts_sha256': digest(run / f'{stem}_contrasts.tsv'),
    'source_regions_sha256': digest(run / f'{stem}_region_grouping.tsv'),
    'promotion_job': sys.argv[2],
    'release_state': 'User adopted corrected Figure 2 estimates on 2026-10-13; COLOC r2 not adopted by this decision',
    'analysis_source': 'figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816',
    'derivative_reason': 'Restore trailing-zero sentinel positions; updated fractions, region counts and interval-only captions',
}
with (run / f'{stem}_manifest.tsv').open('w') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(('field', 'value'))
    writer.writerows(fields.items())
for stem, width, height in (
    ('Fig2D_pip_architecture_by_trait_directness', 1.85, 2.10),
    ('Fig2_companion_noncoding_annotation_gap', 5.50, 2.70),
):
    info = (run / f'{stem}_pdfinfo.txt').read_text()
    assert re.search(r'^Pages:\s+1\s*$', info, re.M), info
    size = re.search(r'^Page size:\s+([\d.]+) x ([\d.]+) pts', info, re.M)
    assert size and abs(float(size[1]) - 72 * width) < 0.01 and abs(float(size[2]) - 72 * height) < 0.01
    fonts = (run / f'{stem}_fonts.txt').read_text()
    assert 'Type 3' not in fonts and 'Helvetica' in fonts, fonts
    artwork = (run / f'{stem}_artwork_text.txt').read_text()
    assert not re.search(r'\b(?:p|q)\s*[<=>]|\bsignificant\b', artwork, re.I), artwork
print('TWO_PANEL_CHECK_PASS: one page each, exact dimensions, Helvetica, no Type 3 or diagnostic significance labels')
PY

"${python_bin}" -B "${project_root}/scripts/figures/validate_fig2d_pip_architecture.py" \
    --panel-dir "${run_root}" > "${run_root}/fig2d_checks.log" 2>&1

"${python_bin}" -B - "${run_root}" <<'PY'
import csv
import hashlib
from pathlib import Path
import sys

run = Path(sys.argv[1])
with (run / 'input_checksums_before.tsv').open() as handle:
    before = list(csv.DictReader(handle, delimiter='\t'))
after = [{'path': row['path'], 'sha256': hashlib.sha256(Path(row['path']).read_bytes()).hexdigest()}
         for row in before]
with (run / 'input_checksums_after.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=('path', 'sha256'), delimiter='\t')
    writer.writeheader()
    writer.writerows(after)
assert before == after, 'Rendering inputs changed during execution'
with (run / 'output_checksums.tsv').open('w') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(('path', 'sha256'))
    for path in sorted(run.iterdir()):
        if path.suffix in {'.pdf', '.tsv', '.txt'} and path.name != 'output_checksums.tsv':
            writer.writerow((path.name, hashlib.sha256(path.read_bytes()).hexdigest()))
(run / 'COMPLETE.txt').write_text('Two individual corrected-coordinate panels rendered; source, PDF and unchanged-input checks passed.\n')
print(f'INPUT_CHECK_PASS: {len(before)} unchanged input/snapshot records')
PY
printf '[fig2] COMPLETE run=%s finished=%s\n' "${run_root}" "$(date -u)"
