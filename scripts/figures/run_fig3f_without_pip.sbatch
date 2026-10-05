#!/usr/bin/env bash
#SBATCH --job-name=fig3f_no_pip
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/fig3f_no_pip_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/logs/fig3f_no_pip_%j.err
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
export PATH=/gpfs/commons/home/jameslee/.local/bin:$PATH
figure_python=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
"$figure_python" scripts/figures/render_fig3f_without_pip.py
figure_out=figures/candidates/fig3f-no-pip-20261005-v1
pdftotext "$figure_out/fig3f_stage_deg_genetics_matrix.pdf" "$figure_out/panel.txt"
pdftoppm -png -r 150 -singlefile "$figure_out/fig3f_stage_deg_genetics_matrix.pdf" "$figure_out/proof"
"$figure_python" - "$figure_out" <<'PY'
from pathlib import Path
import sys,hashlib
from pypdf import PdfReader
from pypdf.generic import ContentStream
out=Path(sys.argv[1]);old=Path('figures/candidates/fig3f-corrected-20261005-v3/source_tables')
for p in old.glob('*.tsv'):
    assert hashlib.sha256(p.read_bytes()).digest()==hashlib.sha256((out/'source_tables'/p.name).read_bytes()).digest()
reader=PdfReader(out/'fig3f_stage_deg_genetics_matrix.pdf');assert len(reader.pages)==1
text=(out/'panel.txt').read_text()
assert 'PIP' not in text and all(t in text for t in ['61 genes','F1-F0','F2-F1','F3-F2','F4-F3','PP.H4'])
fonts=reader.pages[0]['/Resources']['/Font']
for ref in fonts.values():
    font=ref.get_object();assert 'Helvetica' in str(font['/BaseFont'])
    desc=font['/DescendantFonts'][0].get_object()['/FontDescriptor'].get_object()
    assert '/FontFile2' in desc
sizes=[float(args[-1]) for args,op in ContentStream(reader.pages[0].get_contents(),reader).operations if op==b'Tf']
assert sizes and all(abs(s-6)<1e-5 for s in sizes)
(out/'checks.txt').write_text('PASS: no PIP label or column; frozen source tables byte-identical; 61 genes; four stage columns and PP.H4; one page; embedded Helvetica 6 pt.\n')
print((out/'checks.txt').read_text())
PY
