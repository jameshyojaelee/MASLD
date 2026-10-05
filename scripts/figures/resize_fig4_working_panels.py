"""Compact Figure 4D–H using their existing display tables and plot layers.

Run with the rnaseq environment on a compute node. Outputs are new candidates;
the caller installs the reviewed PDFs separately. No scientific models are fit.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from textwrap import dedent

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "figures/main/fig4_singlecell_programs/source_tables/current_candidate"
ML = ROOT / "scripts/analysis/histology_anchored_continuum/molecular_layers"
OUT = Path(sys.argv[1]).resolve()
OUT.mkdir(exist_ok=False)
PANELS = OUT / "panels"
PANELS.mkdir()


def section(code, start, stop):
    return code[code.index(start):code.index(stop)]


def rpath(path):
    return json.dumps(str(path))


# KEY MESSAGE: The same lineage, tissue-state and continuum effects remain
# readable when these individual panels are placed in the manuscript figure.
d = (ROOT / "scripts/figures/render_fig4de_geometry_review.R").read_text()
f = (ML / "64_render_five_cohort_context.R").read_text()
g = (ML / "69_render_refined_molecular_systems_family.R").read_text()
contract = json.loads((ML / "00_contract.json").read_text())
r = '''suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel)})
grDevices::pdf.options(useDingbats=FALSE)
'''
r += f"source({rpath(ROOT / 'scripts/figures/publication_theme.R')})\n"
r += section(d, "theme_panel <- function()", "set.seed(")
r += f"forest <- fread({rpath(SRC / 'fig4d_same_atlas_lineage_specificity.tsv')})\n"
r += '''forest[, lineage_label := factor(comparison_lineage,
  levels=rev(c('Fibroblasts','Cholangiocytes','Endothelial cells','Macrophages','T cells')))]
hero <- unique(forest[, .(program_name=program_key, program_color, program_shape)])
stopifnot(nrow(forest)==10L, sum(forest$target)==2L,
  setequal(hero$program_name, forest$program_key))
'''
r += section(d, "p_d <- ggplot", "bulk <- fread")
r += "stopifnot(all(is.finite(ggplot_build(p_d)$data[[5]]$shape)))\n"
r += f"ggsave({rpath(PANELS / 'fig4d_same_atlas_lineage_specificity.pdf')}, p_d, width=3.3, height=1.85, device=cairo_pdf)\n"
r += "all_cohorts <- c('GSE126848','GSE130970','GSE135251','GSE162694','GSE213621')\n"
r += section(f, 'font_family <-', 'save_panel <- function')
r += section(f, 'cohort_scales <- function', 'make_stage_plot <- function')
r += f"focal_windows <- fread({rpath(SRC / 'fig4f_minimal_20260825/source_tables/focal_hotspot_five_cohort_windows.tsv')})\n"
r += '''figure4f <- make_window_plot(focal_windows, 'focal_program',
  c('Ductular injury','Stromal ECM'), 'Standardized program score',
  show_error=FALSE, legend_labels=cohort_display_plain,
  show_participant_size=FALSE, caption_text=NULL, x_label='Continuum percentile') +
  theme(legend.margin=margin(0,0,0,0), legend.spacing.y=grid::unit(0,'mm')) +
  guides(linetype='none', color=guide_legend(nrow=2, byrow=TRUE))
'''
r += f"ggsave({rpath(PANELS / 'fig4f_continuum_program_trajectories.pdf')}, figure4f, width=3.9, height=2.15, device=cairo_pdf)\n"
r += section(g, 'theme_systems <- function()', 'save_panel <- function')
# Only the system-map label helper is needed; no window ordering is recomputed.
r += '''system_short_label <- function(system_id) sprintf('S%02d', as.integer(sub('system_', '', system_id)))
compact_numeric_labels <- function(values) formatC(values, format='f', digits=1)
'''
r += f"contract <- list(seed={contract['seed']})\n"
r += f"nodes <- fread({rpath(SRC / 'molecular_systems_atlas_20260824/overlays/all_node_stage_continuum_overlays.tsv.gz')})\n"
r += f"communities <- fread({rpath(SRC / 'molecular_systems_atlas_20260824/geometry/community_registry.tsv')})\n"
r += '''system_nodes <- nodes[component_role=='main_system_component']
stopifnot(nrow(system_nodes)==5439L, nrow(communities)==43L)
map_labels <- copy(communities)[order(-n_features)][1:15]
map_labels[, map_label := system_short_label(community_id)]
'''
r += section(g, 'map_label_layer <- function()', 'topology <- ggplot')
r += section(g, 'effect_limits <- function', 'stage_map_path <-')
r = r.replace('save_panel(plot, filename, 4.75, 3.78)',
              'ggsave(file.path(panel_dir, filename), plot, width=3.55, height=3.05, device=cairo_pdf)')
r += f"panel_dir <- {rpath(PANELS)}\n"
r += '''render_effect_map('continuum_beta', 'continuum_supported',
  'Stage-adjusted continuum beta\n(score SD / continuum SD)',
  'Stage-adjusted continuum effects\non the same fixed atlas',
  'fig4g_molecular_systems_continuum_effect_map.pdf')
'''
rfile = OUT / 'render_compact.R'
rfile.write_text(r)
subprocess.run([str(Path(sys.executable).with_name('Rscript')), '--vanilla', str(rfile)], check=True)

sys.path.insert(0, str(ROOT / 'scripts/analysis/fig34_corrections'))
from common import FOCAL, read, plot_setup
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.text import Text

plt = plot_setup()
plt.rcParams.update({'font.family': 'Helvetica', 'pdf.use14corefonts': True,
                     'font.weight': 'normal', 'axes.titleweight': 'normal',
                     'legend.frameon': False, 'axes.unicode_minus': False})
original_save = Figure.savefig


def save(self, *args, **kwargs):
    for ax in self.axes:
        ax.set_title('')
    for text in self.findobj(Text):
        text.set_text(text.get_text().replace('\u2212', '-').replace('\u2013', '-'))
        text.set_fontfamily('Helvetica')
        text.set_fontsize(6)
        text.set_fontweight('normal')
    self.tight_layout(pad=.5)
    return original_save(self, *args, bbox_inches='tight', pad_inches=.04, **kwargs)


Figure.savefig = save
code = (ROOT / 'scripts/analysis/fig34_corrections/05_figures_transport.py').read_text()
e = dedent(section(code, '    # Main 4E', '    # Preserve the actual'))
h = dedent(section(code, '    # Preserve the actual', '    # Cohort and omission'))
# Source estimates and maxT outlines are exactly those in the working panels.
e = e.replace("write(p,out/'fig4e_source.tsv')", "assert len(p)==8")
e = e.replace('figsize=(3.0,2.0)', 'figsize=(2.5,1.8)')
e = e.replace("Cross-sectional contrast versus F0", "Cross-sectional contrast\nversus F0")
# Keep literal newlines inside the Python plotting string escaped.
e = e.replace("contrast\nversus", "contrast\\nversus")
h = h.replace("write(a,out/'fig4h_source.tsv')", "assert len(a)==43 and a.supported.sum()==18")
h = h.replace('figsize=(3.4,2.8)', 'figsize=(2.95,2.45)')
old_labels = section(h, 'offsets=', 'ax.axhline(0')
h = h.replace(old_labels, '''label_positions = {
    'S01':(.286,.705), 'S02':(.04,.46), 'S03':(.036,-.137),
    'S04':(.085,.54), 'S05':(.12,.64), 'S06':(-.10,-.08),
    'S07':(.245,.345), 'S08':(.247,.106), 'S09':(.075,.312),
    'S10':(.148,.265), 'S11':(.285,.24), 'S12':(-.038,.195),
    'S13':(.238,.575), 'S14':(-.035,.29), 'S15':(.28,.47)}
for r in a[a.community_id.isin(labels.community_id)].itertuples():
    ax.annotate(r.figure_label, (r.stage_beta,r.continuum_beta),
                xytext=label_positions[r.figure_label], textcoords='data',
                fontsize=6, arrowprops=dict(arrowstyle='-',lw=.35,color='#777777')
                if r.figure_label in ['S02','S03','S04','S05','S07','S09','S10','S12','S13','S14'] else None)
''')
h = h.replace('Stage coefficient (score SD / cohort stage-code SD)',
              'Stage coefficient\\n(score SD / cohort stage-code SD)')
h = h.replace("fig4h_stage_continuum.pdf", "fig4h_molecular_systems_stage_vs_continuum.pdf")
out = PANELS
exec(compile(e, 'existing_4e_display', 'exec'), globals())
exec(compile(h, 'existing_4h_display', 'exec'), globals())

with (OUT / 'environment.txt').open('w') as handle:
    subprocess.run([str(Path(sys.executable).with_name('Rscript')), '-e', 'sessionInfo()'],
                   check=True, stdout=handle)
    handle.write(f'\nPython {sys.version}\nNumPy {np.__version__}\nPandas {pd.__version__}\n')
inputs = [ROOT / 'scripts/figures/render_fig4de_geometry_review.R',
          ML / '64_render_five_cohort_context.R', ML / '69_render_refined_molecular_systems_family.R',
          ROOT / 'scripts/analysis/fig34_corrections/05_figures_transport.py',
          SRC / 'fig4d_same_atlas_lineage_specificity.tsv',
          SRC / 'fig4f_minimal_20260825/source_tables/focal_hotspot_five_cohort_windows.tsv',
          SRC / 'molecular_systems_atlas_20260824/overlays/all_node_stage_continuum_overlays.tsv.gz',
          SRC / 'molecular_systems_atlas_20260824/geometry/community_registry.tsv',
          ROOT / 'figures/candidates/fig4e-corrected-controls-20260925T175215Z/inputs/fig4e_bulk_tissue_state_transport.tsv',
          SRC / 'molecular_systems_atlas_20260824/overlays/community_stage_continuum_effects.tsv',
          SRC / 'molecular_systems_atlas_20260824/overlays/rendered_label_manifest.tsv',
          SRC / 'molecular_systems_no_delta_20260825/source_tables/system_inference_strip.tsv']
(OUT / 'input_sha256.tsv').write_text('path\tsha256\n' + ''.join(
    f'{p}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n' for p in inputs))
print(OUT)
