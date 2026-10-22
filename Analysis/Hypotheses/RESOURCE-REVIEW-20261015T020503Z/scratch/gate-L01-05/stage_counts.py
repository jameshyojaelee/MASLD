import pandas as pd, numpy as np
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/figures/candidates/bulk-release-corrected-controls-20261008T162308Z/analysis/stage_extensions/'
S='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/gate-L01-05/'
s=pd.read_csv(B+'stage_extension_all_gene_results.tsv',sep='\t',usecols=['gene_id_base','gene_name','logFC','SE','P.Value','FDR','AveExpr','axis','contrast'])
print(s.groupby(['axis','contrast']).size())
f=s[(s.axis=='fibrosis')]
fc=f.pivot_table(index='gene_id_base',columns='contrast',values='logFC')
print('genes F3_vs_F0', fc['F3_vs_F0'].notna().sum(), 'dup gene_id_base:', s[(s.contrast=='F3_vs_F0')&(s.axis=='fibrosis')].gene_id_base.duplicated().sum())
# BH-significant counts per contrast vs Fig 3E (1,338/4,266/6,190/7,043)
for c in ['F1_vs_F0','F2_vs_F0','F3_vs_F0','F4_vs_F0']:
    d=f[f.contrast==c]
    print(c,'FDR<0.05',(d.FDR<0.05).sum(),'FDR<0.05&|lfc|>0.5',((d.FDR<0.05)&(d.logFC.abs()>0.5)).sum())
x=pd.read_csv(S+'classes_rederived.tsv',sep='\t')
f3=f[f.contrast=='F3_vs_F0'].set_index('gene_id_base')
x['in']=x.ensembl.isin(f3.index)
print(x.groupby('cls')['in'].agg(['sum','size']))
x=x[x['in']].copy()
x['mhc']=x.locus.isin(['L090','L091','L092'])
print('MHC genes in table:', x.mhc.sum(), x[x.mhc].cls.value_counts().to_dict())
lc=x.groupby('locus').cls.agg(lambda s: '|'.join(sorted(set(s))))
print('pure-locus counts (incl MHC):'); print(lc[~lc.str.contains(r'\|')].value_counts())
xm=x[~x.mhc]; lcm=xm.groupby('locus').cls.agg(lambda s: '|'.join(sorted(set(s))))
print('pure-locus counts (excl MHC):'); print(lcm[~lcm.str.contains(r'\|')].value_counts())
pure=lcm[~lcm.str.contains(r'\|')]
xp=xm[xm.locus.isin(pure.index)].copy()
xp['logFC']=xp.ensembl.map(f3.logFC); xp['SE']=xp.ensembl.map(f3.SE); xp['AveExpr']=xp.ensembl.map(f3.AveExpr)
gpl=xp.groupby('locus').agg(cls=('cls','first'),ngene=('ensembl','size'),lfc=('logFC','mean'),se2=('SE',lambda v:(v**2).sum()/len(v)**2),ave=('AveExpr','mean'))
print('genes per pure locus by class:'); print(gpl.groupby('cls').ngene.describe())
print('GGT_only vs ALTAST_only loci:', (gpl.cls=='GGT_only').sum(), (gpl.cls=='ALTAST_only').sum())
sub=gpl[gpl.cls.isin(['GGT_only','ALTAST_only'])]
sd=sub.lfc.std(ddof=1)
print('pooled (class-blind) SD of locus-mean F3vsF0 logFC over 2-class loci: %.4f'%sd)
print('mean sampling var share: mean(se^2 of locus mean) %.5f vs total var %.5f'%(sub.se2.mean(), sd**2))
print('AveExpr of locus means: by class (class-level only, no outcome):'); print(sub.groupby('cls').ave.describe())
print('MAD-based SD %.4f'%(1.4826*np.median(np.abs(sub.lfc-sub.lfc.median()))))
sub.to_csv(S+'pure_locus_table_noclassdiff.tsv',sep='\t')
# genome-wide SD of gene logFC
print('genome-wide gene-level F3vsF0 logFC SD %.4f'%f3.logFC.std())
