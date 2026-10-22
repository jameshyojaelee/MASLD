"""Candidate Fig 3F (F3ef, review round 3): released 54-gene selection and colors unchanged;
each released BH mark is drawn filled when it also has BH q<0.05 under joint library-complexity
adjustment (log2 library size, Gini, TMM norm factor; verify2 job 4211) AND with disease-only F0
(f0-reference-arms r8 arm B), and open otherwise.
Adapted from scripts/figures/render_fig3f_without_pip.py (copied as render_fig3f_without_pip.orig.py).

KEY MESSAGE: 26 of the 59 stage marks hold under both checks; 29 of 30 F0-to-F1 marks do not.
Repair round: key for open marks corrected to 'q<0.05, not under both checks'; written to F3ef/repair/.
"""
from pathlib import Path
import hashlib, json, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.font_manager import fontManager
from matplotlib.lines import Line2D
import pandas as pd, numpy as np

root = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
source = root/'figures/candidates/fig3f-no-pip-774-20261012T011128Z/source_tables'
r1 = root/'figures/candidates/review-r3-fixes-20261015T112148Z/F3ef'
out = r1/'repair'
matrix = pd.read_csv(source/'fig3f_stage_genetics_matrix.tsv', sep='\t')
lane = pd.read_csv(source/'fig3f_exact_pair_lanes.tsv', sep='\t').sort_values('row_order')
ordering = pd.read_csv(source/'fig3f_row_order.tsv', sep='\t').sort_values('row_order')
sens = pd.read_csv(r1/'fig3f_mark_sensitivity.tsv', sep='\t')
transitions = ['F0_to_F1','F1_to_F2','F2_to_F3','F3_to_F4']
n = len(ordering); assert n == 54 and len(matrix) == 216
ordered = ordering.ensembl.tolist()
values = matrix.pivot(index='ensembl',columns='contrast',values='logFC').loc[ordered,transitions].to_numpy()
qvalues = matrix.pivot(index='ensembl',columns='contrast',values='FDR').loc[ordered,transitions].to_numpy()
sens['ensembl'] = sens.gene_id_versioned.str.split('.').str[0]
sens['both'] = (sens.q_v4211_joint < .05) & (sens.q_arm_B_disease_only_F0 < .05)
both = sens.pivot(index='ensembl',columns='contrast',values='both').loc[ordered,transitions].to_numpy().astype(bool)
fontManager.addfont('/gpfs/commons/home/jameslee/.local/share/fonts/Helvetica.ttf')
plt.rcParams.update({'font.family':'Helvetica','font.size':6,'pdf.use14corefonts':False,'pdf.fonttype':42,
                     'axes.labelsize':6,'xtick.labelsize':6,'ytick.labelsize':6,'legend.fontsize':6,
                     'axes.unicode_minus':False,'font.weight':'normal','axes.labelweight':'normal'})
fig = plt.figure(figsize=(3.6,n*.105+.95))
grid = fig.add_gridspec(2,2,height_ratios=[1,.035],width_ratios=[4,.55],
                        left=.22,right=.97,top=.97,bottom=.14,wspace=.22,hspace=.12)
axes = [fig.add_subplot(grid[0,k]) for k in range(2)]
lim = float(np.quantile(np.abs(values),.98))
im = axes[0].imshow(values,cmap=LinearSegmentedColormap.from_list('stage',['#1565C0','white','#C9265E']),
                    vmin=-lim,vmax=lim,aspect='auto',interpolation='none')
sig = qvalues < .05
yy,xx = np.where(sig & both)
axes[0].scatter(xx,yy,s=5,c='black',linewidths=0)
yy,xx = np.where(sig & ~both)
axes[0].scatter(xx,yy,s=5,facecolors='none',edgecolors='black',linewidths=.4)
axes[0].set_xticks(range(4),['F1-F0','F2-F1','F3-F2','F4-F3'])
axes[0].set_yticks(range(n),ordering.row_label)
axes[0].set_ylabel(f'{n} genes')
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
handles = [Line2D([],[],marker='o',ls='',ms=2.2,mfc='black',mec='black',mew=0,label=f'q<0.05, also under both checks ({int((sig&both).sum())})'),
           Line2D([],[],marker='o',ls='',ms=2.2,mfc='none',mec='black',mew=.4,label=f'q<0.05, not under both checks ({int((sig&~both).sum())})')]
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,0),ncol=2,frameon=False,handletextpad=.3,columnspacing=1)
pdf = out/'fig3f_stage_deg_genetics_matrix_mark_sensitivity.pdf'
fig.savefig(pdf,bbox_inches='tight',pad_inches=.04)
plt.close(fig)
src = out/'fig3f_source_tables'; src.mkdir(exist_ok=True)
for p in source.glob('*.tsv'): shutil.copy(p, src/p.name)
summ = {'genes':n,'effects':4*n,'marks':int(sig.sum()),'marks_filled':int((sig&both).sum()),'marks_open':int((sig&~both).sum()),
        'filled_by_contrast':dict(zip(transitions,(sig&both).sum(0).tolist())),'marks_by_contrast':dict(zip(transitions,sig.sum(0).tolist())),
        'genes_with_filled_mark':int((sig&both).any(1).sum()),'models_refit':False}
(out/'fig3f_mark_sensitivity_summary.json').write_text(json.dumps(summ,indent=2)+'\n'); print(summ)
