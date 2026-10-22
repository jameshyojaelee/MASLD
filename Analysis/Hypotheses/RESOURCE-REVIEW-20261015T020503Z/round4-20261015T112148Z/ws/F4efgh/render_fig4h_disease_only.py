#!/usr/bin/env python3
"""Candidate Figure 4H on the disease-only stage-complete roster (109 GSE162694 + 258 GSE213621).

KEY MESSAGE: Among disease units only, system fibrosis-stage and stage-adjusted continuum
coefficients still align (most systems share direction), less tightly than when the 68
GSE213621 controls entered the stage trend as stage 0.

Same marks, size, colors and typography as the installed panel
(figures/candidates/fig3-fig4-layout-774-20261013T135430Z/producers/render_fig4h_single_outline.py).
Differences: x and y come from fig4h_disease_only/run_20261015b/system_stage_continuum_by_roster.tsv
(roster disease_only, x = stage, y = cont_fixed_projection); no outline is drawn, because the
size-matched competitive null exists only for the 435-unit roster that includes the controls;
x-axis label names the fibrosis-stage coefficient. No statistic is fitted here.
Usage: render_fig4h_disease_only.py OUT_DIR
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
WS = ROOT / 'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/ws/F4efgh'
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=False)
font_manager.fontManager.addfont('/gpfs/commons/home/jameslee/.local/share/fonts/Helvetica.ttf')
plt.rcParams.update({'font.family': 'Helvetica', 'font.size': 6, 'pdf.fonttype': 42,
                     'pdf.use14corefonts': False, 'axes.unicode_minus': False,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': .5, 'xtick.major.width': .5, 'ytick.major.width': .5})

sysb = pd.read_csv(WS / 'fig4h_disease_only/run_20261015b/system_stage_continuum_by_roster.tsv', sep='\t')
src = ROOT / 'figures/main/fig4_singlecell_programs/source_tables'
inst = pd.read_csv(src / 'corrected_774/fig4h_single_outline/fig4h_source.tsv', sep='\t')
labels = pd.read_csv(src / 'current_candidate/molecular_systems_atlas_20260824/overlays/rendered_label_manifest.tsv', sep='\t')
labels = labels[(labels.figure_id == 'continuum_effect_map') &
                (labels.selection_rule == '15_largest_systems_fixed_before_outcome_overlay')]
assert len(labels) == 15
d = sysb[sysb.roster == 'disease_only'][['community_id', 'stage', 'cont_fixed_projection']].rename(
    columns={'stage': 'stage_beta', 'cont_fixed_projection': 'continuum_beta'})
pub = sysb[sysb.roster == 'as_published'].set_index('community_id')
chk = inst.set_index('community_id').loc[pub.index]
assert np.allclose(pub.stage, chk.stage_beta, atol=1e-12) and np.allclose(pub.cont_fixed_projection, chk.continuum_beta, atol=1e-12)
a = d.merge(inst[['community_id', 'figure_label']], on='community_id', validate='one_to_one')
assert len(a) == 43
a['stage_beta_installed_435'] = a.community_id.map(chk.stage_beta)
a['continuum_beta_installed_435'] = a.community_id.map(chk.continuum_beta)
a['direction_pair'] = np.where((a.stage_beta > 0) & (a.continuum_beta > 0), 'Both positive',
                      np.where((a.stage_beta < 0) & (a.continuum_beta < 0), 'Both negative', 'Opposite signs'))
a['fixed_outcome_blind_label'] = a.community_id.isin(labels.community_id)
rho = spearmanr(a.stage_beta, a.continuum_beta)[0]
same = int((np.sign(a.stage_beta) == np.sign(a.continuum_beta)).sum())
a['roster'] = 'disease_only_367_units_109_GSE162694_258_GSE213621'
a['x_unit'] = 'outcome_score_per_cohort_SD_of_recorded_fibrosis_stage_disease_units'
a['y_unit'] = 'outcome_score_SD_per_continuum_SD_adjusted_for_recorded_stage_category_and_sex_disease_units'
a.to_csv(out / 'fig4h_disease_only_source.tsv', sep='\t', index=False)

cols = {'Both positive': '#C9265E', 'Both negative': '#2166AC', 'Opposite signs': '#9E9E9E'}
fig, ax = plt.subplots(figsize=(215 / 72, 176 / 72))
for lab, col in cols.items():
    g = a[a.direction_pair == lab]
    ax.scatter(g.stage_beta, g.continuum_beta, c=col, s=15, linewidths=0, zorder=3)
ax.axhline(0, lw=.4, color='#9E9E9E', zorder=1)
ax.axvline(0, lw=.4, color='#9E9E9E', zorder=1)
offs = {'S15': (-4, 8), 'S07': (8, 1), 'S04': (6, -7), 'S14': (-16, 6)}
for r in a[a.fixed_outcome_blind_label].itertuples():
    dx, dy = offs.get(r.figure_label, (3, 2))
    ax.annotate(r.figure_label, (r.stage_beta, r.continuum_beta), xytext=(dx, dy),
                textcoords='offset points', fontsize=6,
                arrowprops=dict(arrowstyle='-', lw=.35, color='#777777', shrinkA=0, shrinkB=2)
                if r.figure_label in offs else None)
ax.set_xlabel('Fibrosis-stage coefficient, disease units\n(score SD / cohort stage-code SD)', fontsize=6)
ax.set_ylabel('Stage-adjusted continuum coefficient\n(score SD / continuum SD)', fontsize=6)
ax.tick_params(labelsize=6, length=2, pad=1)
ax.legend(handles=[Line2D([], [], marker='o', ls='', ms=3, color=c, label=l) for l, c in cols.items()],
          frameon=False, loc='upper left', fontsize=6, handletextpad=.2, borderaxespad=.1, labelspacing=.25)
ax.text(.98, .02, 'Spearman rho = %.2f;\n%d/43 same direction' % (rho, same), transform=ax.transAxes,
        ha='right', va='bottom', fontsize=6, linespacing=1.2)
fig.tight_layout(pad=.5)
pdf = out / 'fig4h_molecular_systems_stage_vs_continuum_disease_only.pdf'
fig.savefig(pdf)
plt.close(fig)
print('rho', rho, 'same', same, 'pdf', pdf)
