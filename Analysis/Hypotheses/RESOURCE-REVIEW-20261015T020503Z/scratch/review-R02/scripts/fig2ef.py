import pandas as pd, numpy as np
B='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
S=B+'Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/scratch/review-R02/'
m=pd.read_pickle(S+'r2_main.pkl')
reg=pd.read_csv(B+'GWAS/finemapping/config/gwas_registry.tsv',sep='\t')
anc=dict(zip(reg.study_name,reg.ancestry))
m['anc']=m.gwas_name.map(anc)
print('registry vs master ancestry mismatch:',(m.anc!=m.ancestry).sum())
m=m[m.gene.notna()&(m.gene.astype(str).str.strip()!='')]
sig=m[m['PP.H4.susie']>0.5]
bg=sig.groupby('gene').anc.apply(set)
def binof(a):
    if a=={'EUR'}: return 'EUR only'
    if 'EUR' in a: return 'shared'
    return 'nonEUR'
gb=bg.apply(binof)
print(gb.value_counts(), len(gb))
non=gb[gb=='nonEUR'].index
# locus-level check: for each non-EUR-unique gene, its best SuSiE row; nearest EUR SuSiE-positive row top_snp within 1Mb
sig=sig.copy()
sig['c']=sig.top_snp.str.split(':').str[0]; sig['p']=sig.top_snp.str.split(':').str[1].astype(int)
eur=sig[sig.anc=='EUR']
abfeur=m[(m.anc=='EUR')&(m['PP.H4.abf']>0.5)].copy()
abfeur['c']=abfeur.top_snp.str.split(':').str[0]; abfeur['p']=abfeur.top_snp.str.split(':').str[1].astype(int)
rows=[]
for g in non:
    r=sig[sig.gene==g].sort_values('PP.H4.susie',ascending=False).iloc[0]
    e=eur[(eur.c==r.c)&((eur.p-r.p).abs()<=1_000_000)]
    a=abfeur[(abfeur.c==r.c)&((abfeur.p-r.p).abs()<=1_000_000)]
    same=e[e.top_snp==r.top_snp]
    rows.append(dict(gene=g,gwas=r.gwas_name,snp=r.top_snp,pp4=round(r['PP.H4.susie'],3),
      eur_susie_genes_1Mb=','.join(sorted(e.gene.unique()))[:80], n_eur_susie_1Mb=e.gene.nunique(),
      eur_same_snp=','.join(sorted(same.gene.unique())), n_eur_abf_1Mb=a.gene.nunique()))
t=pd.DataFrame(rows)
pd.set_option('display.width',250)
print(t.to_string())
print('non-EUR-unique genes with a EUR SuSiE-positive gene within 1 Mb:',(t.n_eur_susie_1Mb>0).sum(),'of',len(t))
print('... with EUR SuSiE or ABF positive within 1Mb:',((t.n_eur_susie_1Mb>0)|(t.n_eur_abf_1Mb>0)).sum())
t.to_csv(S+'fig2e_locus_check.tsv',sep='\t',index=False)
# 2F
d=m.copy()
bs=d[d['PP.H4.susie'].notna()].sort_values('PP.H4.susie',ascending=False).groupby(['gene','anc']).head(1)[['gene','anc','PP.H4.susie','gwas_name','top_snp']]
ba=d[d['PP.H4.abf'].notna()].sort_values('PP.H4.abf',ascending=False).groupby(['gene','anc']).head(1)[['gene','anc','PP.H4.abf','gwas_name','top_snp']]
c=bs.merge(ba,on=['gene','anc'],how='outer',suffixes=('_s','_a'))
c['sup']=(c['PP.H4.susie']>0.5)|(c['PP.H4.abf']>0.5)
c['multi']=c['PP.H4.susie']>0.5
r=c.groupby('gene').agg(n_sup=('sup','sum'),n_multi=('multi','sum'))
g13=r[r.n_sup>=4].sort_values(['n_multi','n_sup'],ascending=False)
print(len(g13), list(g13.index), (r.n_sup>=5).sum())
cc=c[c.gene.isin(g13.index)&c.sup].copy()
cc['lead']=np.where(cc.multi,cc.top_snp_s,cc.top_snp_a); cc['study']=np.where(cc.multi,cc.gwas_name_s,cc.gwas_name_a)
for g in g13.index:
    x=cc[cc.gene==g]
    print(g, '|', '; '.join(f"{a}:{s}:{l}" for a,s,l in zip(x.anc,x.study,x.lead)))
# how many ancestry-supported cells rest on single-signal ABF only
print('cells in 13 genes:',len(cc),'ABF-only cells:',(~cc.multi).sum())
