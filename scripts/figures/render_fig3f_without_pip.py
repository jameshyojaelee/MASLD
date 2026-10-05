"""Remove the PIP display using the frozen, already-reviewed 3F tables.

KEY MESSAGE: Display the same stage coefficients and exact-pair PP.H4.
"""
from pathlib import Path
import hashlib
import json
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.font_manager import fontManager
import pandas as pd
import numpy as np

root = Path(__file__).resolve().parents[2]
source = root/'figures/candidates/fig3f-corrected-20261005-v3/source_tables'
out = root/'figures/candidates/fig3f-no-pip-20261005-v1'
out.mkdir(exist_ok=False)
shutil.copytree(source, out/'source_tables')
matrix = pd.read_csv(source/'fig3f_stage_genetics_matrix.tsv', sep='\t')
lane = pd.read_csv(source/'fig3f_exact_pair_lanes.tsv', sep='\t').sort_values('row_order')
ordering = pd.read_csv(source/'fig3f_row_order.tsv', sep='\t').sort_values('row_order')
transitions = ['F0_to_F1','F1_to_F2','F2_to_F3','F3_to_F4']
assert len(ordering)==61 and len(matrix)==244 and len(lane)==61
ordered = ordering.ensembl.tolist()
values = matrix.pivot(index='ensembl',columns='contrast',values='logFC').loc[ordered,transitions].to_numpy()
qvalues = matrix.pivot(index='ensembl',columns='contrast',values='FDR').loc[ordered,transitions].to_numpy()
fontManager.addfont('/gpfs/commons/home/jameslee/.local/share/fonts/Helvetica.ttf')
plt.rcParams.update({'font.family':'Helvetica','font.size':6,'pdf.use14corefonts':False,'pdf.fonttype':42,
                     'axes.labelsize':6,'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':6,
                     'axes.unicode_minus':False,'font.weight':'normal','axes.labelweight':'normal'})
fig = plt.figure(figsize=(3.6,61*.105+.75))
grid = fig.add_gridspec(2,2,height_ratios=[1,.035],width_ratios=[4,.55],
                        left=.22,right=.97,top=.97,bottom=.12,wspace=.22,hspace=.12)
axes = [fig.add_subplot(grid[0,k]) for k in range(2)]
lim = float(np.quantile(np.abs(values),.98))
im = axes[0].imshow(values,cmap=LinearSegmentedColormap.from_list('stage',['#1565C0','white','#C9265E']),
                    vmin=-lim,vmax=lim,aspect='auto',interpolation='none')
yy,xx = np.where(qvalues < .05)
axes[0].scatter(xx,yy,s=1.6,c='black',linewidths=0)
axes[0].set_xticks(range(4),['F1-F0','F2-F1','F3-F2','F4-F3'])
axes[0].set_yticks(range(61),ordering.row_label)
axes[0].set_ylabel('61 genes')
pp = axes[1].imshow(lane[['PP.H4']].to_numpy(),cmap=LinearSegmentedColormap.from_list('pph4',['white','#ED5565']),
                    vmin=0,vmax=1,aspect='auto',interpolation='none')
axes[1].set_xticks([0],['PP.H4']);axes[1].set_yticks([])
for ax in axes:
    ax.tick_params(length=0,pad=3)
    for spine in ax.spines.values(): spine.set_visible(False)
cb = fig.colorbar(im,cax=fig.add_subplot(grid[1,0]),orientation='horizontal',ticks=[-lim,0,lim])
cb.set_label('log2FC');cb.set_ticklabels([f'{-lim:.1f}','0',f'{lim:.1f}'])
cb.outline.set_visible(False);cb.ax.tick_params(length=2,width=.4)
cb2 = fig.colorbar(pp,cax=fig.add_subplot(grid[1,1]),orientation='horizontal',ticks=[0,1])
cb2.outline.set_visible(False);cb2.ax.tick_params(length=2,width=.4)
pdf = out/'fig3f_stage_deg_genetics_matrix.pdf'
fig.savefig(pdf,bbox_inches='tight',pad_inches=.04);plt.close(fig)
rows=[{'role':'frozen_input','path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(source.glob('*.tsv'))]
rows += [{'role':role,'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
         for role,p in [('renderer',Path(__file__)),('output_pdf',pdf)]]
pd.DataFrame(rows).to_csv(out/'manifest.tsv',sep='\t',index=False)
(out/'summary.json').write_text(json.dumps({'genes':61,'effects':244,'change':'PIP column and crosses removed',
                                         'models_refit':False,'source_tables_unchanged':True},indent=2)+'\n')
print(out)
