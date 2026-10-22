"""Figure 4E on the 774-unit roster: same marks, size and typography as the
working panel (producers/render_DE.py of fig4-assembly-review-20261005-v1).
Reads the fitted program_axis_map.tsv only; no model is refitted here."""
from pathlib import Path
import json
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
from pypdf import PdfReader

# Copied from figures/candidates/fig4e-program-axis-774-20261012T010830Z/producers/render_fig4e_774.py.
# Changes: renders only the 774-unit three-cohort arm into the review-r3 candidate directory,
# and adds the truncated (ad hoc) Knapp-Hartung call counts from the review recompute to panel_checks.
SRC_ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/fig4e-program-axis-774-20261012T010830Z')
OUT = Path(sys.argv[1])
WORKING = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/main/'
               'fig4_singlecell_programs/panels/fig4e_bulk_tissue_state_transport.pdf')
ARMS = {
    'new774_three_cohorts': dict(
        panel=OUT / 'panels',
        ylim=(0, 2.4),
        state='774_unit_roster_one_record_per_genetic_identity_cluster; three cohorts with fibrosis and NAS'),
}
font_manager.fontManager.addfont('/gpfs/commons/home/jameslee/.local/share/fonts/Helvetica.ttf')
plt.rcParams.update({'font.family': 'Helvetica', 'font.size': 6, 'axes.labelsize': 6,
                     'axes.titlesize': 6, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                     'legend.fontsize': 6, 'pdf.fonttype': 42, 'pdf.use14corefonts': False,
                     'axes.unicode_minus': False, 'axes.spines.top': False,
                     'axes.spines.right': False})
XLIM = (-.3, .45)
checks = {}
for arm, cfg in ARMS.items():
    axis = pd.read_csv(SRC_ROOT / 'analysis' / arm / 'program_axis_map.tsv', sep='\t')
    assert len(axis) == 117 and axis.feature_id.nunique() == 117
    eligible = axis[axis.testable].copy()
    assert (eligible.weight_coverage >= .8 - 1e-12).all()
    assert (axis.loc[~axis.testable, 'weight_coverage'] < .8).all()
    # Independent BH over the pooled fibrosis and NAS family.
    ps = np.concatenate([eligible['p_value_meta_' + a].to_numpy() for a in ['fibrosis', 'nas']])
    order = np.argsort(ps)
    q = np.empty(len(ps))
    q[order] = np.minimum(1, np.minimum.accumulate(
        (ps[order] * len(ps) / np.arange(1, len(ps) + 1))[::-1])[::-1])
    stored = np.concatenate([eligible['q_value_' + a].to_numpy() for a in ['fibrosis', 'nas']])
    np.testing.assert_allclose(q, stored, rtol=1e-10, atol=1e-12)
    cfg['panel'].mkdir(parents=True, exist_ok=False)
    source = cfg['panel'].parent / 'source_tables'
    source.mkdir(exist_ok=False)
    axis.to_csv(source / 'fig4e_all117_bulk_axis_map.tsv', sep='\t', index=False)
    long = []
    for a in ['fibrosis', 'nas']:
        for row in eligible.itertuples():
            long.append(dict(
                program_uid=row.feature_id, axis=a, beta=getattr(row, 'beta_meta_' + a),
                ci_low=getattr(row, 'ci_lower_' + a), ci_high=getattr(row, 'ci_upper_' + a),
                q=getattr(row, 'q_value_' + a), n_cohorts=getattr(row, 'n_cohorts_' + a),
                effect_unit='cohort-standardized program score SD / ' +
                ('recorded stage-code unit' if a == 'fibrosis' else 'NAS point'),
                comparison='partial linear association; sex and other histology axis adjusted',
                evidence_state=cfg['state'],
                interval_display='retained in source table; volcano encodes effect and BH q'))
    long = pd.DataFrame(long)
    long.to_csv(source / 'fig4e_plotted_effects.tsv', sep='\t', index=False)
    assert long.beta.between(*XLIM).all(), 'an effect falls outside the x axis'
    assert (-np.log10(long.q)).between(*cfg['ylim']).all(), 'a q value falls outside the y axis'

    fig, axes = plt.subplots(1, 2, figsize=(180 / 72, 108 / 72), sharey=True)
    fig.subplots_adjust(left=.19, right=.97, bottom=.27, top=.88, wspace=.24)
    for ax, a, title, unit in zip(axes, ['fibrosis', 'nas'], ['Fibrosis: linear', 'NAS: linear'],
                                  ['SD / stage-code\nunit', 'SD / NAS point']):
        beta = eligible['beta_meta_' + a].to_numpy()
        qq = eligible['q_value_' + a].to_numpy()
        colors = np.where(qq < .05, np.where(beta > 0, '#C9265E', '#1565C0'), '#9E9E9E')
        idx = np.argsort(qq < .05)
        ax.scatter(beta[idx], -np.log10(qq[idx]), s=4, c=colors[idx], linewidths=0)
        ax.axhline(-np.log10(.05), color='#707070', ls=(0, (2, 2)), lw=.4)
        ax.axvline(0, color='#D0D0D0', lw=.3)
        ax.set(xlim=XLIM, ylim=cfg['ylim'], xticks=[-.2, 0, .4], yticks=[0, 1, 2])
        ax.set_title(title, pad=2)
        ax.set_xlabel(unit, labelpad=1)
        ax.tick_params(length=2, pad=1, width=.4)
        for spine in ax.spines.values():
            spine.set_linewidth(.4)
    axes[0].set_ylabel('-log10 BH q', labelpad=1)
    pdf = cfg['panel'] / 'fig4e_bulk_tissue_state_transport.pdf'
    fig.savefig(pdf)
    plt.close(fig)

    reader = PdfReader(pdf)
    page = reader.pages[0]
    working = PdfReader(WORKING).pages[0]
    fonts = []
    for ref in page['/Resources'].get('/Font', {}).values():
        f = ref.get_object()
        desc = f
        if '/DescendantFonts' in f:
            desc = f['/DescendantFonts'][0].get_object()
        fd = desc.get('/FontDescriptor')
        fd = fd.get_object() if fd is not None else {}
        fonts.append(dict(base=str(f.get('/BaseFont')), subtype=str(f.get('/Subtype')),
                          descendant_subtype=str(desc.get('/Subtype')),
                          embedded=any(k in fd for k in ('/FontFile', '/FontFile2', '/FontFile3'))))
    text = subprocess.run(['pdftotext', '-layout', str(pdf), '-'], check=True,
                          capture_output=True, text=True).stdout
    subprocess.run(['pdftoppm', '-r', '300', '-singlefile', '-png', str(pdf),
                    str(pdf.with_suffix(''))], check=True)
    subprocess.run(['gs', '-q', '-dNOPAUSE', '-dBATCH', '-sDEVICE=nullpage', str(pdf)], check=True)
    checks[arm] = dict(
        pdf=str(pdf), pages=len(reader.pages),
        page_size_pt=[float(page.mediabox.width), float(page.mediabox.height)],
        working_page_size_pt=[float(working.mediabox.width), float(working.mediabox.height)],
        page_size_equals_working=(float(page.mediabox.width) == float(working.mediabox.width) and
                                  float(page.mediabox.height) == float(working.mediabox.height)),
        fonts=fonts,
        all_fonts_helvetica_embedded=bool(fonts) and all('Helvetica' in f['base'] and f['embedded'] for f in fonts),
        no_type3=all('Type3' not in (f['subtype'] + f['descendant_subtype']) for f in fonts),
        pdftotext=text.split(),
        n_registry=117, n_eligible=int(len(eligible)), bh_family=int(2 * len(eligible)),
        n_points_per_facet=int(len(eligible)),
        calls_fibrosis=int((eligible.q_value_fibrosis < .05).sum()),
        calls_fibrosis_positive=int(((eligible.q_value_fibrosis < .05) & (eligible.beta_meta_fibrosis > 0)).sum()),
        calls_fibrosis_negative=int(((eligible.q_value_fibrosis < .05) & (eligible.beta_meta_fibrosis < 0)).sum()),
        calls_nas=int((eligible.q_value_nas < .05).sum()),
        max_neglog10_q=float((-np.log10(long.q)).max()),
        beta_range=[float(long.beta.min()), float(long.beta.max())],
        ylim=list(cfg['ylim']), xlim=list(XLIM))
    assert checks[arm]['pages'] == 1 and checks[arm]['page_size_equals_working']
    assert checks[arm]['all_fonts_helvetica_embedded'] and checks[arm]['no_type3']
(OUT / 'panel_checks.json').write_text(json.dumps(checks, indent=2) + '\n')
print(json.dumps(checks, indent=2))
