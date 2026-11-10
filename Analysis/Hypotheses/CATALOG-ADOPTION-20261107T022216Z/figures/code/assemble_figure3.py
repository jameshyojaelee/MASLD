"""KEY MESSAGE: discovery and held-out stage effects use the current individual panels.

Candidate only. Panels remain vector, at native scale; this explicitly requested
review assembly does not overwrite the Illustrator source or canonical assembly.
"""
from collections import Counter
from pathlib import Path
import csv
import hashlib
import json
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pypdf import PdfReader, PdfWriter

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
OUT = ROOT / 'Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/figures/review_assembly'
PANELS = ROOT / 'figures/main/fig3_bulk_transcriptomics/panels'
font_manager.fontManager.addfont('/gpfs/commons/home/jameslee/.local/share/fonts/Helvetica.ttf')
plt.rcParams.update({'font.family': 'Helvetica', 'font.size': 6, 'pdf.fonttype': 42})
W, H = 636, 987
# x, top, display label, source file; no scaling.
SLOTS = [
    (8, 26, 'A', 'fig3a_cohort_metadata_matrix.pdf'),
    (232, 26, 'B', 'fig3b_nas_fib_grid.pdf'),
    (8, 212, 'C', 'fig3c_cohort_alluvial.pdf'),
    (8, 456, 'D', 'fig3d_pca_fibrosis_gradient.pdf'),
    (8, 680, 'E', 'fig3e_stage_remodeling.pdf'),
    (340, 585, 'E, held-out replication', 'fig3e_steatosite_replication_strip.pdf'),
    (370, 26, 'F', 'fig3f_stage_deg_genetics_matrix.pdf'),
    (370, 753, 'G', 'fig3g_genetic_state_complementarity.pdf'),
]
target = OUT / 'Figure3.pdf'
assert not target.exists(), target
writer = PdfWriter()
page = writer.add_blank_page(width=W, height=H)
records, boxes = [], []
for x, top, label, name in SLOTS:
    path = (ROOT / 'Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/figures/panels' / name) if name.startswith('fig3g') else PANELS / name
    reader = PdfReader(path)
    assert len(reader.pages) == 1
    p = reader.pages[0]
    w, h = float(p.mediabox.width), float(p.mediabox.height)
    y = H - top - h
    assert x >= 0 and y >= 0 and x + w <= W and y + h <= H
    for bx, by, bw, bh, blabel in boxes:
        assert x+w <= bx or bx+bw <= x or y+h <= by or by+bh <= y, (label, blabel)
    boxes.append((x, y, w, h, label))
    page.merge_translated_page(p, x, y, expand=False)
    records.append(dict(label=label, source=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        x=x, y=y, width=w, height=h, scale=1))
fig = plt.figure(figsize=(W/72, H/72), facecolor='none')
ax = fig.add_axes([0, 0, 1, 1])
ax.set(xlim=(0, W), ylim=(0, H)); ax.axis('off')
for x, top, label, _ in SLOTS:
    ax.text(x, H-top+5, label, fontsize=6, fontweight='normal', va='bottom')
overlay = OUT / 'Figure3_letters.pdf'
fig.savefig(overlay, transparent=True); plt.close(fig)
page.merge_page(PdfReader(overlay).pages[0])
page.compress_content_streams()
with target.open('wb') as f:
    writer.write(f)
def text(path):
    return subprocess.run(['pdftotext', '-layout', str(path), '-'], check=True, text=True, capture_output=True).stdout
actual = text(target)
expected = Counter()
for row in records:
    expected += Counter(''.join(text(Path(row['source'])).split()))
missing = expected - Counter(''.join(actual.split()))
assert not missing, missing
assert 'SteatoSITE' in actual and '1,267' in actual and '6 cohorts' in actual
(OUT / 'Figure3_layout.txt').write_text(actual)
with (OUT / 'Figure3_panel_sources.tsv').open('w') as f:
    z = csv.DictWriter(f, fieldnames=records[0], delimiter='\t'); z.writeheader(); z.writerows(records)
(OUT / 'Figure3_checks.json').write_text(json.dumps(dict(pages=1, width=W, height=H,
    transform='translation_only', panel_characters_missing=0, current_3A_present=True), indent=2)+'\n')
print(target)
