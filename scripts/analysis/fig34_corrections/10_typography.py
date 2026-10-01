"""Render frozen candidate estimates with Helvetica, uniform 6 pt, and no titles."""
import sys, importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
from common import ROOT, FOCAL, PROGRAM, read, write, plot_setup
from scipy import stats

base = Path(sys.argv[1]).resolve()
dest = base / 'typography_v2'
version = 2
while dest.exists():
    version += 1
    dest = base / f'typography_v{version}'
dest.mkdir(exist_ok=False)
for name in ['source_tables_v2', 'lineage_matched', 'continuum']:
    (dest / name).symlink_to(base / name, target_is_directory=True)
out = dest / 'figures_v4'
out.mkdir()
plt = plot_setup()
plt.rcParams.update({'font.family': 'Helvetica', 'pdf.use14corefonts': True,
                     'font.weight': 'normal', 'axes.titleweight': 'normal',
                     'legend.title_fontsize': 6, 'legend.frameon': False,
                     'axes.unicode_minus': False})
from matplotlib.figure import Figure
from matplotlib.text import Text
original_save = Figure.savefig
def save(self, *args, **kwargs):
    # Cohort/program identity belongs in axis labels, not plot headings.
    for ax in self.axes:
        title = ax.get_title()
        if title and not title.startswith('Historical'):
            ax.set_xlabel(title + '\n' + ax.get_xlabel())
        ax.set_title('')
        ax.set_facecolor('white')
    if self._suptitle is not None:
        self._suptitle.set_text('')
    for text in self.findobj(Text):
        text.set_text(text.get_text().replace('\u2212', '-').replace('\u2013', '-'))
        text.set_fontfamily('Helvetica')
        text.set_fontsize(6)
        text.set_fontweight('normal')
    self.tight_layout(pad=.5)
    kwargs.setdefault('bbox_inches', 'tight')
    kwargs.setdefault('pad_inches', .04)
    return original_save(self, *args, **kwargs)
Figure.savefig = save

# KEY MESSAGE: Retain the original effects and memberships; correct display only.
code = (Path(__file__).parent / '05_figures_transport.py').read_text()
start = code.index('    # Main 4E')
stop = code.index('    # Previously evaluated cohorts')
from textwrap import dedent
exec(dedent(code[start:stop]), globals())
gs = ROOT/'figures/candidates/fig3g-coloc-r2-20260924T061535Z/source_tables'
scatter = read(gs/'fig3g_raw_treat_interface.tsv')
fig, ax = plt.subplots(figsize=(3.4, 2.8))
for state, col in [('Neither','#9E9E9E'),('Transcriptomics sig.','#C9265E'),
                   ('Genetics sig.','#1565C0'),('Both','#6A3D9A')]:
    g = scatter[scatter.state.eq(state)]
    ax.scatter(g.logFC,g.susie_pp4,s=3,c=col,alpha=.7,
               label=f'{state} ({len(g)})',linewidths=0)
ax.set_xlabel('Disease-state log2FC'); ax.set_ylabel('SuSiE PP.H4')
ax.legend(frameon=False)
fig.savefig(out/'fig3g_genetics_state_interface.pdf'); plt.close(fig)
pn = read(base/'protein/protein_pair_covariance.tsv')
fig, axs = plt.subplots(1,3,figsize=(5.5,2.7),sharex=True)
for ax,(pair,g) in zip(axs,pn.groupby('pair',sort=False)):
    ax.errorbar(g.rho,np.arange(len(g)),xerr=[g.rho-g.low,g.high-g.rho],
                fmt='o',ms=3,lw=.6,color='#1565C0')
    ax.set_yticks(range(len(g)),g.model.str.replace('_',' '))
    ax.set_xlabel(pair+'\nPartial r (95% CI)')
    ax.axvline(0,color='#9E9E9E',lw=.5)
fig.savefig(out/'S3_protein_covariance.pdf'); plt.close(fig)

# Reuse the existing display tables, never refit scientific models.
finish = (Path(__file__).parent/'09_panel_finish.py').read_text()
finish = finish.replace('figsize=(7.2,2.05)', 'figsize=(5.5,2.2)')
finish = finish.replace('figsize=(6.8,2.8)', 'figsize=(5.5,3.0)')
finish = finish.replace(";write(test,base/'source_tables_v2/six_hallmark_testability.tsv')", '')
finish = finish.replace("ax.text(.5,.7,'EPITHELIAL MESENCHYMAL TRANSITION',ha='center',fontsize=7);", '')
finish = finish.replace('Untestable: 150/200', 'EMT untestable: 150/200')
# Normalize only the PDF text check; source-category joins retain exact labels.
finish = finish.replace('all(v in txt for v in', "all(v.replace('\u2013', '-') in txt.replace('\u2013', '-') for v in")
namespace = {'__name__': 'typography_finish'}
exec(compile(finish,'09_panel_finish.py','exec'),namespace)
namespace['main'](dest)
(dest/'CAPTION_NOTES.txt').write_text(
    'Titles removed. Bulk effects remain conditional on unresolved donor identity.\n'
    'Protein example is selection-conditioned; same 46 participants.\n'
    'Content-lineage effects use matched correction-available cells and donors.\n'
    'Figure 3G: PP.H4 is unsigned; disease-state log2FC is signed.\n'
    'Only display styling changed; no scientific models refit.\n')
print(dest)
