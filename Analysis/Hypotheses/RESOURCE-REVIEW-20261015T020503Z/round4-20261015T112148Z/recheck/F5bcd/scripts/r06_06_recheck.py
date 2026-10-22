import pandas as pd, numpy as np
from scipy import stats
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
prot=pd.read_csv(B+'Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv',sep='\t')
print(prot.columns.tolist(), len(prot))
deg=pd.read_csv(B+'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/deg_results.csv')
print(deg.columns.tolist()[:12], len(deg))
atlas=pd.read_csv(B+'RNA-seq/results/multi_evidence/multi_evidence_atlas.csv',usecols=['human_symbol','ensembl_id'])
deg['ens']=deg.gene.str.replace(r'\..*$','',regex=True)
m=atlas.rename(columns={'human_symbol':'g','ensembl_id':'ens'}).merge(deg[['ens','symbol']],on='ens')
m=m[m.g==m.symbol][['g','ens']]
cnt=deg[deg.symbol.notna()&(deg.symbol!='')].symbol.value_counts()
uniq=set(cnt[cnt==1].index)
rep=deg[deg.symbol.isin(uniq)&~deg.symbol.isin(m.g)][['symbol','ens']].rename(columns={'symbol':'g'})
m=pd.concat([m,rep])
assert not m.g.duplicated().any() and not m.ens.duplicated().any()
bulk=m.merge(deg[['ens','logFC','padj']],on='ens')
con=prot[['gene','logFC','padj']].rename(columns={'gene':'g','logFC':'plfc','padj':'ppadj'}).merge(bulk,on='g')
con=con[np.isfinite(con.logFC)&np.isfinite(con.plfc)]
print('n full',len(con),'rho full',stats.spearmanr(con.logFC,con.plfc)[0])
rd=lambda f: set(x.strip() for x in open(B+'Analysis/Spatial/results/universe_validation/'+f) if x.strip())
uni=rd('prioritized_universe_FINAL.txt'); gen=rd('universe_genetic.txt'); tx=rd('universe_transcriptomic.txt')
print('uni',len(uni),'gen',len(gen),'tx',len(tx),'uni==gen|tx',uni==(gen|tx))
du=con[con.g.isin(uni)].copy()
print('n uni',len(du),'rho uni',stats.spearmanr(du.logFC,du.plfc)[0],'conc',(np.sign(du.logFC)==np.sign(du.plfc)).mean())
du['tx']=du.g.isin(tx); du['gen']=du.g.isin(gen)
print(du.groupby(['tx','gen']).size())
for lab,s in [('tx arm',du[du.tx]),('gen only',du[~du.tx&du.gen]),('gen arm any',du[du.gen]),('neither',du[~du.tx&~du.gen])]:
    if len(s)>2: print(lab,len(s),round(stats.spearmanr(s.logFC,s.plfc)[0],4))
ps=du[du.ppadj<0.05]
print('protein sig',len(ps))
for lab,s in [('rna<0.05',ps[ps.padj<0.05]),('rna>=0.05',ps[~(ps.padj<0.05)])]:
    a=(np.sign(s.logFC)==np.sign(s.plfc)).sum(); n=len(s)
    print(lab,a,n,a/n,stats.binomtest(a,n,0.5).pvalue)
print('NA rna padj among protein sig', ps.padj.isna().sum())
print('not protein sig', (~(du.ppadj<0.05)).sum())
inst=pd.read_csv(B+'figures/main/fig5_molecular_context/panels/data/fig5b_protein_triage_genes.tsv',sep='\t')
mm=inst.merge(du,left_on='gene',right_on='g')
print('installed rows',len(inst),'matched',len(mm),'max |dLFC|',(mm.bulk_logFC-mm.logFC).abs().max(),(mm.protein_logFC-mm.plfc).abs().max())
