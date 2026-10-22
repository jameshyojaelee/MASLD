#!/usr/bin/env python3
"""Candidate S4 focal continuum forest on the 774-unit disease-only evaluation roster.

KEY MESSAGE: Stromal ECM and Ductular injury rise along the fixed projection after stage and
sex adjustment in both evaluation cohorts; adding cell-type composition leaves Ductular injury
near its unadjusted size, while the Stromal ECM estimate falls by an amount that depends on
which cell-type fraction set is used.

Replaces the installed panel, which prints the 844-record roster (GSE213621 n=293).
Inputs (no model is fitted here):
  basis_774/continuum/focal_effects.tsv (rebuilt fractions) and
  sensitivity_774_scores_earlier_fractions/continuum/focal_effects.tsv (earlier fractions),
  both in figures/candidates/continuum-composition-774-20261012T012137Z; rows with
  axis == fixed_projection and arm == disease_only. Same drawing as
  scripts/analysis/fig34_corrections/09_panel_finish.py (continuum forest), plus the
  earlier-fraction row, without the suptitle.
Usage: render_s4_focal_continuum_774.py OUT_DIR
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

BASE = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/continuum-composition-774-20261012T012137Z')
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=False)
font_manager.fontManager.addfont('/gpfs/commons/home/jameslee/.local/share/fonts/Helvetica.ttf')
plt.rcParams.update({'font.family': 'Helvetica', 'font.size': 6, 'pdf.fonttype': 42,
                     'pdf.use14corefonts': False, 'axes.unicode_minus': False,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .5})
FOCAL = {'hotspot_hepatocytes_f05c535ae5bbc0b9': 'Stromal ECM', 'hotspot_hepatocytes_48f39dd4d817a10e': 'Ductular injury'}
rows = []
for src, tag in [('basis_774', 'rebuilt'), ('sensitivity_774_scores_earlier_fractions', 'earlier')]:
    c = pd.read_csv(BASE / src / 'continuum/focal_effects.tsv', sep='\t')
    c = c[(c.axis == 'fixed_projection') & (c.arm == 'disease_only') & c.feature_id.isin(FOCAL)].copy()
    c['fraction_set'] = np.where(c.composition, tag, 'none')
    c['source'] = str(BASE / src / 'continuum/focal_effects.tsv')
    rows.append(c if tag == 'rebuilt' else c[c.composition])
t = pd.concat(rows, ignore_index=True)
order = [('GSE162694', 'none'), ('GSE162694', 'rebuilt'), ('GSE162694', 'earlier'),
         ('GSE213621', 'none'), ('GSE213621', 'rebuilt'), ('GSE213621', 'earlier')]
lab = {'none': 'Stage/sex adjusted', 'rebuilt': 'Composition adjusted,\nrebuilt fractions',
       'earlier': 'Composition adjusted,\nearlier fractions'}
colr = {'none': '#2166AC', 'rebuilt': '#C9265E', 'earlier': '#E89AB4'}
assert len(t) == 12
assert set(t[t.dataset == 'GSE213621'].n) == {258} and set(t[t.dataset == 'GSE162694'].n) == {109}
keep = ['feature_id', 'dataset', 'fraction_set', 'n', 'beta', 'hc3_se', 'hc3_low', 'hc3_high', 'hc3_p', 'hc3_q', 'family_size', 'source']
t['program'] = t.feature_id.map(FOCAL)
t[['program'] + keep].to_csv(out / 'S4_focal_continuum_cohort_effects_774_source.tsv', sep='\t', index=False)

fig, axs = plt.subplots(1, 2, figsize=(396 / 72, 230 / 72))
for ax, (uid, title) in zip(axs, FOCAL.items()):
    g = t[t.feature_id == uid].set_index(['dataset', 'fraction_set']).loc[order[::-1]]
    for j, ((ds, fs), r) in enumerate(g.iterrows()):
        ax.errorbar(r.beta, j, xerr=[[r.beta - r.hc3_low], [r.hc3_high - r.beta]], fmt='o', ms=3, lw=.6,
                    color=colr[fs], capsize=0)
    ax.set_yticks(range(len(g)), [f'{ds} ({lab[fs]}), n={int(r.n)}' if False else f'{ds}, n={int(r.n)}\n{lab[fs]}'
                                  for (ds, fs), r in g.iterrows()], fontsize=6)
    ax.axvline(0, lw=.5, color='#9E9E9E')
    ax.axhline(2.5, lw=.3, color='#D0D0D0')
    ax.set_title(title, fontsize=6, pad=2)
    ax.set_xlabel('Program SD / continuum SD\n(HC3 95% CI)', fontsize=6)
    ax.tick_params(labelsize=6, length=2, pad=1)
fig.tight_layout(pad=.5)
pdf = out / 'S4_focal_continuum_cohort_effects_774.pdf'
fig.savefig(pdf)
plt.close(fig)
print(t[['program', 'dataset', 'fraction_set', 'n', 'beta', 'hc3_low', 'hc3_high', 'hc3_q']].to_string())
