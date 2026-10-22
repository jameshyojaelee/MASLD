"""NG17 gate counts: coverage, positive control, in-data effect estimate, MDE."""
import numpy as np, pandas as pd, re
from scipy.stats import norm
rng = np.random.default_rng(20260923)
R='/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/'
P=R+'figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/'
g=pd.read_parquet(P+'passport_gene_index.parquet'); g['ens']=g.ensembl_id.astype(str).str.split('.').str[0]
print('catalog rows',len(g),'unique ens',g.ens.nunique())
print('route counts',g.next_experiment_rule_id.value_counts(dropna=False).to_dict())
print('class counts',g.primary_evidence_class.value_counts(dropna=False).to_dict())
st=g[g.next_experiment_rule_id=='EXP_STRESS_PERTURB_V1']
print('stress rows',len(st),'unique ens',st.ens.nunique()); print('stress classes',st.primary_evidence_class.value_counts().to_dict())
for c in g.columns:
    if g[c].dtype==object and g[c].nunique()<12: pass
print('gene_index cols',list(g.columns))
stress=set(st.ens)
gn=pd.read_csv(R+'data/external/gnomad_constraint/gnomad.v4.1.constraint_metrics.tsv',sep='\t',
  usecols=['gene','gene_id','mane_select','lof.oe_ci.upper','lof.pLI','syn.oe','cds_length'])
gn=gn[gn.mane_select.astype(str).str.lower().isin(['true','1'])].copy()
gn['ens']=gn.gene_id.astype(str).str.split('.').str[0]
gn=gn.dropna(subset=['lof.oe_ci.upper']).drop_duplicates('ens')
gn['dec']=pd.qcut(gn['lof.oe_ci.upper'].rank(method='first'),10,labels=False)
gn['bottom']=(gn.dec==0).astype(int)
print('MANE genes with LOEUF',len(gn),'bottom decile LOEUF cut',gn.loc[gn.bottom==1,'lof.oe_ci.upper'].max())
hpa=pd.read_csv(R+'data/external/hpa/proteinatlas.tsv',sep='\t',usecols=['Gene','Ensembl','Secretome location','Blood concentration - Conc. blood IM [pg/L]','Blood concentration - Conc. blood MS [pg/L]'])
hpa=hpa.drop_duplicates('Ensembl').set_index('Ensembl')
# DepMap
dm=pd.read_csv(R+'Analysis/downstream_analysis/essentiality/CRISPRGeneEffect.csv',index_col=0)
print('DepMap lines',dm.shape[0],'genes',dm.shape[1])
entrez={int(re.search(r'\((\d+)\)',c).group(1)):c for c in dm.columns}
sym={c.split(' (')[0]:c for c in dm.columns}
# common essential approx: gene effect < -0.5 in >=90% of lines (Chronos); and DepMap-style 90th percentile rank
frac=(dm< -0.5).mean(0)
ce_thr=set(frac[frac>=0.9].index)
rk=dm.rank(axis=1,pct=True)  # low effect = low rank
ce_rank=set(rk.columns[(rk.quantile(0.9,axis=0))<=0.1]) # gene in bottom decile of each line at 90th percentile line
print('common-essential approx: threshold-rule',len(ce_thr),' rank-rule',len(ce_rank))
# map ens->DepMap via gnomAD symbol (current)
gn['dm']=gn.gene.map(sym)
print('MANE genes mapped to DepMap',gn.dm.notna().sum())
# coverage
def cov(S,name):
    S=set(S); x=gn[gn.ens.isin(S)]
    h=hpa.reindex(list(S))
    print(f'{name}: n={len(S)} LOEUF={len(x)} pLI={x["lof.pLI"].notna().sum()} DepMap={x.dm.notna().sum()} HPA={h.Gene.notna().sum()} '
          f'secretome_annot={h["Secretome location"].notna().sum()} bloodIM={h["Blood concentration - Conc. blood IM [pg/L]"].notna().sum()} bloodMS={h["Blood concentration - Conc. blood MS [pg/L]"].notna().sum()}')
cov(stress,'stress')
# controls: unresolved
co=pd.read_csv(R+'RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/primary_source_controls/coefficient_table.tsv.gz',sep='\t')
co['ens']=co.gene.str.split('.').str[0]; co=co.drop_duplicates('ens').set_index('ens')
gn['ave']=gn.ens.map(co.AveExpr)
for cls in g.primary_evidence_class.dropna().unique():
    pass
unres_classes=[c for c in g.primary_evidence_class.dropna().unique() if 'unresolv' in str(c).lower() or 'insufficient' in str(c).lower() or 'no_' in str(c).lower()]
print('candidate unresolved classes',unres_classes)
unres=set(g[g.primary_evidence_class.isin(unres_classes)].ens)-stress
d=gn[gn.ave.notna()].copy()
d['grp']=np.where(d.ens.isin(stress),1,np.where(d.ens.isin(unres),0,-1))
d=d[d.grp>=0].copy()
d['aved']=pd.qcut(d.ave.rank(method='first'),10,labels=False)
d['lend']=pd.qcut(d.cds_length.rank(method='first'),10,labels=False)
print('analysable stress',(d.grp==1).sum(),'unresolved pool',(d.grp==0).sum())
d['ce']=d.dm.isin(ce_thr).astype(int)
for v in ['bottom','ce']:
    print(f'crude {v}: stress {d.loc[d.grp==1,v].mean():.4f} unresolved {d.loc[d.grp==0,v].mean():.4f}')
# stratified (matched) estimate: reweight controls to case strata
d['s']=d.aved.astype(str)+'_'+d.lend.astype(str)
cs=d[d.grp==1].s.value_counts()
ctrl=d[d.grp==0]; cc=ctrl.s.value_counts()
ok=cs.index[cs.index.isin(cc.index)]
print('case strata',len(cs),'with controls',len(ok),'cases in strata with >=3x controls',int(cs[[s for s in ok if cc[s]>=3*cs[s]]].sum()))
for v in ['bottom','ce']:
    exp=sum(cs[s]*ctrl.loc[ctrl.s==s,v].mean() for s in ok)/cs[ok].sum()
    obs=d[(d.grp==1)&d.s.isin(ok)][v].mean()
    print(f'matched {v}: stress {obs:.4f} expected {exp:.4f} diff {obs-exp:+.4f}')
# matched permutation, 2000 draws (gate only)
cases=d[(d.grp==1)&d.s.isin(ok)]
groups={s:ctrl.loc[ctrl.s==s,'bottom'].values for s in ok}
obs=cases.bottom.sum()
null=np.array([sum(rng.choice(groups[s],cs[s],replace=True).sum() for s in ok) for _ in range(2000)])
print('perm bottom: obs',obs,'null mean',null.mean(),'two-sided p',min(1,2*min((1+(null>=obs).sum())/2001,(1+(null<=obs).sum())/2001)))
# MDE: two-proportion, n1 cases vs 3n1 controls, alpha .05 two-sided, power .8
def mde(n1,p0,k=3):
    n0=k*n1
    for dlt in np.arange(0.001,0.5,0.0005):
        p1=p0+dlt; pb=(n1*p1+n0*p0)/(n1+n0)
        se0=np.sqrt(pb*(1-pb)*(1/n1+1/n0)); se1=np.sqrt(p1*(1-p1)/n1+p0*(1-p0)/n0)
        if (dlt-norm.ppf(.975)*se0)/se1>=norm.ppf(.8): return dlt
n1=len(cases)
for v,p0 in [('bottom',ctrl.bottom.mean()),('ce',ctrl.ce.mean())]:
    print(f'MDE {v}: n1={n1} p0={p0:.4f} MDE(+)={mde(n1,p0):.4f} (alpha .05); holm-level alpha .0167:',end=' ')
    n0=3*n1; dl=np.arange(0.001,0.5,0.0005)
    z=norm.ppf(1-0.0167/2)
    for x in dl:
        p1=p0+x; pb=(n1*p1+n0*p0)/(n1+n0)
        if (x-z*np.sqrt(pb*(1-pb)*(1/n1+1/n0)))/np.sqrt(p1*(1-p1)/n1+p0*(1-p0)/n0)>=norm.ppf(.8): print(round(x,4)); break
# positive control CEG2 / negative NEG1
ceg=pd.read_csv(R+'data/public_r4/Hart2017_CEG2/CEGv2.txt',sep='\t'); neg=pd.read_csv(R+'data/public_r4/Hart2017_CEG2/NEGv1.txt',sep='\t')
def map_list(df):
    cols=[entrez.get(int(e)) for e in df.ENTREZ_ID if pd.notna(e) and str(e).isdigit() or isinstance(e,(int,np.integer))]
    return set(c for c in cols if c)
print('NEG cols',list(neg.columns))
cegc=map_list(ceg)
negc=map_list(neg.rename(columns={neg.columns[0]:'GENE'})) if 'ENTREZ_ID' in neg.columns else set(sym[s] for s in neg.iloc[:,0] if s in sym)
print('CEG2 n',len(ceg),'mapped DepMap',len(cegc),'| NEG n',len(neg),'mapped',len(negc))
gn['ceg']=gn.dm.isin(cegc); gn['neg']=gn.dm.isin(negc)
for lab,m in [('CEG2',gn.ceg),('NEG1',gn.neg)]:
    x=gn[m]; print(f'{lab} with LOEUF {len(x)} bottom-decile frac {x.bottom.mean():.3f} (background 0.100); median LOEUF {x["lof.oe_ci.upper"].median():.3f} vs all {gn["lof.oe_ci.upper"].median():.3f}; common-ess(thr) {x.dm.isin(ce_thr).mean():.3f}; syn.oe median {x["syn.oe"].median():.3f} vs {gn["syn.oe"].median():.3f}')
print('stress genes common-essential (thr)',d.loc[d.grp==1,'ce'].sum())
