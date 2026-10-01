"""Typography, control colors and explicit untestable panel; no scientific fits."""
import sys,shutil,subprocess,hashlib
import numpy as np
import pandas as pd
from common import ROOT,ML,FOCAL,read,write,plot_setup

def main(base):
    out=base/'panels';out.mkdir(exist_ok=False);qa=base/'verification_final';qa.mkdir(exist_ok=False);plt=plot_setup()
    for p in (base/'figures_v4').glob('*.pdf'):
        if p.name.startswith(('S4G','S4N','S3_hallmark')) or p.name=='S4_content_lineage_disease_effects.pdf':continue
        shutil.copy2(p,out/p.name)
    for prefix in ['S4G','S4N','S3_hallmark_stage']:
        tab=read(base/f'source_tables_v2/{prefix}_stage_display.tsv')
        for uid,g in tab.groupby('feature_id',sort=False):
            name=FOCAL.get(uid,uid.removeprefix('HALLMARK_').replace('_',' '))
            fig,axs=plt.subplots(1,4,figsize=(7.2,2.05),sharey=True)
            for ax,co in zip(axs,['GSE130970','GSE135251','GSE162694','GSE213621']):
                order=['Control','F0–F1','F2','F3–F4'] if co=='GSE213621' else ['F0','F1','F2','F3','F4']
                z=g[g.dataset.eq(co)].set_index('fibrosis_group_reported').reindex(order)
                for j,r in enumerate(z.itertuples()):
                    if r.count>0:
                        ax.errorbar(j,r.mean,yerr=1.96*r.std/np.sqrt(r.count),fmt='o',ms=3,lw=.6,color='#9E9E9E' if order[j]=='Control' else '#C9265E')
                    else:ax.text(j,.04,'No score',transform=ax.get_xaxis_transform(),rotation=90,va='bottom',ha='center',fontsize=5,color='#777777')
                ax.set_xticks(range(len(order)),order);ax.set_xlim(-.3,len(order)-.7);ax.set_title(co+('\nGrouped histology' if co=='GSE213621' else '\nExact stage'));ax.axhline(0,lw=.4,color='#9E9E9E')
                if co=='GSE213621':ax.set_facecolor('#F3F3F3')
            axs[0].set_ylabel('Mean stored score');fig.suptitle(name+' | Historical roster; identity conditional',fontsize=6);fig.tight_layout();fig.savefig(out/f'{prefix}_{uid}.pdf');plt.close(fig)
    test=read(ML/'figures/source_tables/s3_six_hallmark_testability.tsv');write(test,base/'source_tables_v2/six_hallmark_testability.tsv')
    missing=test.groupby('set_id').testable.any();missing=missing[~missing].index
    assert list(missing)==['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']
    fig,ax=plt.subplots(figsize=(4.4,1.5));ax.axis('off');ax.text(.5,.7,'EPITHELIAL MESENCHYMAL TRANSITION',ha='center',fontsize=7);ax.text(.5,.4,'Untestable: 150/200 original members observed (75%)\nFixed minimum coverage: 80%; family remains 50 Hallmarks',ha='center',va='center',fontsize=6);fig.tight_layout();fig.savefig(out/'S3_hallmark_stage_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf');plt.close(fig)
    l=read(base/'lineage_matched/lineage_cohort_effects.tsv');l=l[l.minimum_cells.eq(50)&l.annotation_filter.eq('all')&l.contrast.str.startswith('primary')&l.arm.ne('raw_full')]
    c=read(base/'continuum/focal_effects.tsv');c=c[c.axis.eq('fixed_projection')&c.arm.eq('disease_only')]
    for kind,tab,name in [('lineage',l,'S4_content_lineage_disease_effects'),('continuum',c,'S4_focal_continuum_cohort_effects')]:
        fig,axs=plt.subplots(1,2,figsize=(6.8,2.8))
        for ax,(uid,title) in zip(axs,FOCAL.items()):
            g=tab[tab['program_uid' if kind=='lineage' else 'feature_id'].eq(uid)]
            for j,r in enumerate(g.itertuples()):ax.errorbar(r.beta,j,xerr=[[r.beta-r.hc3_low],[r.hc3_high-r.beta]],fmt='o',ms=3,lw=.6,color=('#9E9E9E' if r.arm=='raw' else '#C9265E') if kind=='lineage' else '#2166AC')
            if kind=='lineage':labels=[f'{r.dataset}\n'+('Paired raw' if r.arm=='raw' else 'Corrected')+f', n={r.n}' for r in g.itertuples()];xlabel='MASH − steatosis projected score\n(HC3 95% CI)'
            else:labels=[f'{r.dataset}\n'+('Composition adjusted' if r.composition else 'Stage/sex adjusted')+f', n={r.n}' for r in g.itertuples()];xlabel='Program SD / continuum SD\n(HC3 95% CI)'
            ax.set_yticks(range(len(g)),labels);ax.set_xlabel(xlabel);ax.set_title(title);ax.axvline(0,lw=.5,color='#9E9E9E')
        fig.suptitle('Same correction-available cells and donors' if kind=='lineage' else 'Historical library roster; identity conditional',fontsize=6);fig.tight_layout();fig.savefig(out/(name+'.pdf'));plt.close(fig)
    texts=qa/'pdf_text';texts.mkdir();proofs=qa/'proofs';proofs.mkdir();manifest=[]
    for p in sorted(out.glob('*.pdf')):
        txt=subprocess.check_output(['pdftotext','-raw',str(p),'-'],text=True);assert len(txt.strip())>20
        (texts/(p.stem+'.txt')).write_text(txt)
        if p.stem.startswith(('S4G','S4N','S3_hallmark')) and 'EPITHELIAL_MESENCHYMAL_TRANSITION' not in p.stem:assert all(v in txt for v in ['Control','F0–F1','F3–F4','Grouped histology'])
        subprocess.run(['pdftoppm','-scale-to','1400','-singlefile','-png',str(p),str(proofs/p.stem)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
        manifest.append(dict(panel=p.stem,path=str(p.resolve()),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),candidate_only=True,pdf_text_check=True))
    assert len(manifest)==25
    write(pd.DataFrame(manifest),qa/'figure_manifest.tsv')
    (qa/'AUTOMATED_CHECKS_COMPLETE').write_text('25 PDF text/category checks pass; estimates reused unchanged.\n')

if __name__=='__main__':
    from pathlib import Path
    main(Path(sys.argv[1]))
