#!/usr/bin/env bash
#SBATCH --job-name=fig3f_corrected
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/fig3f_corrected_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/fig3f_corrected_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
export PATH=/gpfs/commons/home/jameslee/.local/bin:$PATH
project_root=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
figure_out="$project_root/figures/candidates/fig3f-corrected-20261005-v3"
python_bin=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$python_bin" scripts/figures/render_fig3f_corrected.py --out "$figure_out"
mkdir "$figure_out/proofs"
pdftotext "$figure_out/fig3f_stage_deg_genetics_matrix.pdf" "$figure_out/proofs/panel.txt"
pdftoppm -png -r 180 -singlefile "$figure_out/fig3f_stage_deg_genetics_matrix.pdf" "$figure_out/proofs/panel"
"$python_bin" - "$figure_out" <<'PY'
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import pdfplumber
from pypdf import PdfReader
from pypdf.generic import ContentStream
out = Path(sys.argv[1])
summary = json.loads((out/'summary.json').read_text())
matrix = pd.read_csv(out/'source_tables/fig3f_stage_genetics_matrix.tsv', sep='\t')
lane = pd.read_csv(out/'source_tables/fig3f_exact_pair_lanes.tsv', sep='\t')
assert len(matrix) == 4*summary['genes'] and len(lane) == summary['genes']
assert matrix.groupby('ensembl').FDR.min().lt(.05).all()
assert lane.eligible.eq(True).all() and lane['PP.H4'].gt(.5).all()
assert lane.pip.isna().all()
# Independently compare every displayed PP.H4 with the completed gene-study summary.
root = Path.cwd()
studies = pd.read_csv(root/'GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv',
                      usecols=['ensembl','gwas_name','PP.H4.susie']).dropna(subset=['ensembl'])
check = lane.merge(studies[['ensembl','gwas_name','PP.H4.susie']], left_on=['ensembl','gwas'], right_on=['ensembl','gwas_name'], validate='one_to_one')
assert len(check) == len(lane) and np.max(np.abs(check['PP.H4']-check['PP.H4.susie'])) < 1e-8
pdf = out/'fig3f_stage_deg_genetics_matrix.pdf'
reader = PdfReader(pdf)
assert len(reader.pages) == 1
fonts = [str(font.get_object()['/BaseFont']) for font in reader.pages[0]['/Resources']['/Font'].values()]
assert set(fonts) == {'/Helvetica'}, fonts
(out/'proofs/fonts.txt').write_text('\n'.join(fonts)+'\n')
sizes = [float(args[-1]) for args, op in ContentStream(reader.pages[0].get_contents(), reader).operations if op == b'Tf']
assert sizes and all(abs(size-6) < 1e-5 for size in sizes), set(sizes)
with pdfplumber.open(pdf) as p:
    page = p.pages[0]
    assert all(c['x0'] >= 0 and c['x1'] <= page.width and c['top'] >= 0 and c['bottom'] <= page.height for c in page.chars)
text = (out/'proofs/panel.txt').read_text()
assert all(x in text for x in ['F1-F0','F2-F1','F3-F2','F4-F3','PP.H4','PIP','log2FC'])
assert 'Historical' not in text and 'Figure 3F' not in text
(out/'checks.txt').write_text('PASS: full BH families, exact-pair and gene-study agreement, missing-PIP state, 4 effects per gene, single-page Helvetica normal 6 pt, required labels, no title, no clipped text.\n')
print((out/'checks.txt').read_text())
PY
